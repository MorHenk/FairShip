from ROOT import TVector3, TRotation, TMatrixD, TVectorD, TDecompLU
import shipunit as u

def calculate_center_from_lot(list_of_tubes):
    '''Calculates the geometric center point of a DtModule from a given list of DriftTube objects
    that are passed to this function as parameter list_of_tubes.
    
    Parameters
    ----------
    list_of_tubes: list
        List of drift tube objects that are part of the module for which the geometric center is to be computed
        
    Returns
    -------
    TVector3
        Vector to the center point of the module
    '''
    x_sum = 0
    y_sum = 0
    z_sum = 0
    for tube in list_of_tubes:
        tube_center = tube.get_center_position()
        x_sum += tube_center[0]
        y_sum += tube_center[1]
        z_sum += tube_center[2]
    
    center_x = x_sum / len(list_of_tubes)
    center_y = y_sum / len(list_of_tubes)
    center_z = z_sum / len(list_of_tubes)
    
    return TVector3(center_x,center_y,center_z)
    
def calculate_rotation_from_lot(list_of_tubes):
    return -1 
    
def parse_det_id(det_id):
    """Parse the detector id to human readable dictionary so that a specific tube can easily be
    identified and addressed out of a bigger detector arrangement
    
    Parameters
    ----------
    det_id: int
        Detector ID that is to be parsed
        
    Returns
    -------
    dict
        dictionary containing the result in human readable form
    """
    result = {}
    str_id = str(det_id)
    last_two = int(str_id[-2:]) #last four digits of the detector ID
    #parse view
    view = "X"
    result['station'] = det_id / 10000000
    if result['station'] == 1:
        if str_id[1] == '1':
            view = "U"
    elif result['station'] == 2:
        if str_id[1] == '0':
            view = "V"
            
    #parse module
    module = "T" + str(result['station'])
    if result['station'] >= 3:
        if last_two <= 12:
            module += "d"
        elif last_two <= 24:
            module += "c"
        elif last_two <= 36:
            module += "b"
        else:
            module += "a"
            
    module += view
    result['module'] = module
    return result

def calculate_center(vec1, vec2):
    '''Calculate the position in the middle between two positions pointed to by vectors.
    
    Parameters
    ----------
    vec1: TVector3
        Vector to the first position
        
    vec2: TVector3
        Vector to the second position
        
    Returns
    -------
    TVector3
        Position of the point in between the two positions passed as arguments
    '''
    vec_between = vec2 - vec1
    half_length = vec_between.Mag() / 2
    vec_between.SetMag(half_length)
    
    return TVector3(vec1 + vec_between)
    
def z_rotation_to_euler_angles(rad_z):
    '''Convert a rotation around the z-axis to Euler angles following the X-convention.
    
    Parameters
    ----------
    rad_z: float
        Radians of rotation around z
        
    Returns
    -------
    float,float,float
        Phi, Theta, Psi following the Euler-X-convention
    '''
    rot = TRotation()
    rot.RotateZ(rad_z)    
    return rot.GetXPhi(), rot.GetXTheta(), rot.GetXPsi()

def distance_to_wire(tube,mom=None,pos=None):
    """Calculates the distance of closest approach for a track and a tube in mm.
    Note: This distance is positive if a valid track was used.
    
    Parameters
    ----------
    tube: DriftTube
        DriftTube object that was hit.
        
    mom: ROOT.TVector3
        Momentum (a.k.a direction) of the track
    
    pos: ROOT.TVector3
        Position on the track
        
    Returns
    -------
    float
        Closest distance between track and wire
    """
    vtop,vbot = tube.wire_end_positions()    
    normal_vector = mom.Cross(vtop-vbot)
    vec_any_two_points = vbot - pos #a vector from any point on one of the straight lines to any point on the other straight
    distance = (abs(vec_any_two_points.Dot(normal_vector)) / normal_vector.Mag())

    return distance

def measurement_vector(tube,mom,pos):
    """ Calculates the vector which is perpendicular to the track AND the sense wire of the passed tube. This is, when the passed track
    is used as seed for a genfit.GBL refit, the prediction for the measurement at that specific sense wire.
    
    
    Parameters
    ----------
    tube: DriftTube
        DriftTube object that was hit
        
    mom: ROOT.TVector3
        Momentum (a.k.a direction) of the track
    
    pos: ROOT.TVector3
        Position on the track 
    
    Returns
    -------
    ROOT.TVector3
        Vector of closest approach in lab system (x,y,z). Note, that this contains no information about the point of closest approach itself.
        Unit of returned vector is mm.
    """
    vtop,vbot = tube.wire_end_positions()
    wire_dir = vtop - vbot
    plane = {'pos' : TVector3(pos - vbot),
             'dir1' : TVector3(mom),
             'dir2' : TVector3(-1 * wire_dir)
        }
      
    const_vector = TVectorD(2)
    coeff_matrix = TMatrixD(2,2)
    
    const_vector[0] = -(plane['pos'].Dot(mom))
    const_vector[1] = -(plane['pos'].Dot(wire_dir))
    
    coeff_matrix[0][0] = plane['dir1'].Dot(mom)
    coeff_matrix[0][1] = plane['dir2'].Dot(mom)
    coeff_matrix[1][0] = plane['dir1'].Dot(wire_dir)
    coeff_matrix[1][1] = plane['dir2'].Dot(wire_dir)
    
    #Solve linear eqs:
    solvable_matrix = TDecompLU(coeff_matrix)
    result = TVectorD(const_vector)
    rc = solvable_matrix.Solve(result)
    if not rc:
        raise Exception("Matrix is not solvable")
    
    #Extrapolate track and wire by the result times their directions
    PCA_on_track = TVector3(pos + result[0] * mom)
    PCA_on_wire = TVector3(vbot + result[1] * wire_dir)
    
    return TVector3(PCA_on_track - PCA_on_wire)
    

def calculate_residuals(track,dtmodules,module_residuals):
    """ Calculates the residuals for a given track and returns these in a dictionary
    grouped to modules of 48 drift tubes.
    
    Parameters
    ----------
    track: 
        genfit Track object for that the residuals should be calculated
    dtmodules:
        dictionary containing the drift tube modules as DtAlignment.DtModule objects
    module_residuals:
        dictionary containing the residuals per module with the module name as keys.
        This is where the result is written to
    """     
    if __debug__:
        # Check for conistency
        for key in dtmodules.keys():
            if not key in module_residuals.keys():
                print("Error: key {} not in residuals dictionary".format(key))
                
    # 1.) Loop over hits in track
    n_points = track.getNumPointsWithMeasurement()
    points = track.getPointsWithMeasurement()
    for i in range(n_points):
        point = points[i]
        raw_measurement = point.getRawMeasurement()
        det_id = raw_measurement.getDetId()
        rt_dist = raw_measurement.getRawHitCoords()[6] #rt distance stored here
        # 2.) Parse detector id to module
        module_id = parse_det_id(det_id)
        module = dtmodules[module_id['module']]
        # 3.) Find correct drift tube in module
        for j in range(len(module.get_tubes())):
            tube = module.get_tubes()[j]
            if tube._ID == det_id:
                break
        tube = module.get_tubes()[j]
    
        # 4.) Read fitted position and momentum from fitted state
        fitted_state = track.getFittedState(i)
        mom = fitted_state.getMom()
        pos = fitted_state.getPos()
        
        # 5.) Calculate distance of track to wire
        dist = distance_to_wire(tube, mom, pos)
        # 6.) Calculate residual and append to correct entry in dictionary
        residual = dist - rt_dist
        module_residuals[module_id['module']].append(residual)
        
def calculate_residuals_lr(track,dtmodules,module_residuals):
    n_points = track.getNumPointsWithMeasurement()
    points = track.getPointsWithMeasurement()
    for i in range(n_points):
        point = points[i]
        raw_measurement = point.getRawMeasurement()
        det_id = raw_measurement.getDetId()
        rt_dist = raw_measurement.getRawHitCoords()[6] #rt distance stored here
        # 2.) Parse detector id to module
        module_id = parse_det_id(det_id)
        module = dtmodules[module_id['module']]
        # 3.) Find correct drift tube in module
        for j in range(len(module.get_tubes())):
            tube = module.get_tubes()[j]
            if tube._ID == det_id:
                break
        tube = module.get_tubes()[j]
    
        # 4.) Read fitted position and momentum from fitted state
        fitted_state = track.getFittedState(i)
        mom = fitted_state.getMom()
        pos = fitted_state.getPos()
        
        # 5.) Calculate distance of track to wire
        dist = measurement_vector(tube, mom, pos)
        # 6.) Calculate residual and append to correct entry in dictionary
        residual = dist.Mag() - rt_dist
        rot = TRotation()
        rot.SetXAxis(dist)
        res_vec = TVector3(residual,0,0)
        backrot = rot * res_vec
        if backrot.X() > 0: #missed left, towards pos x
            module_residuals[module_id['module']]['l'].append(residual)
        else:
            module_residuals[module_id['module']]['r'].append(residual)        
