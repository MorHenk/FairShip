#ifndef CHARMDET_MILLEPEDECALLER_H_
#define CHARMDET_MILLEPEDECALLER_H_

#include "TObject.h"
#include "Mille.h"
//includes for GBL fitter from genfit
#include <vector>
#include "GblPoint.h"
#include "GblTrajectory.h"
#include "MilleBinary.h"
#include "Track.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include <map>
#include "TVector3.h"
#include "TDecompLU.h"
#include "TRotation.h"
#include "TMath.h"
#include <cstdint>


/**
 * A class for wrapping the millepede function call such that it can be called from
 * within a python script
 *
 * @author Stefan Bieschke
 * @date Apr. 9, 2019
 */
class MillepedeCaller//: public TObject
{
public:
	MillepedeCaller(const char *outFileName, bool asBinary = true, bool writeZero = false);
	~MillepedeCaller();

	void call_mille(int n_local_derivatives,
					const float *local_derivatives,
					int n_global_derivatives,
					const float *global_derivatives,
					const int *label,
					float measured_residual,
					float sigma);

	const int* labels() const;
	double perform_GBL_refit(const genfit::Track& track) const;

	ClassDef(MillepedeCaller,3);

private:
	Mille mille;
	gbl::MilleBinary* m_gbl_mille_binary;

	//helper methods
	std::vector<gbl::GblPoint> list_hits(const genfit::Track* track) const;
	TMatrixD* calc_jacobian(const genfit::Track* track, const unsigned int hit_id_1, const unsigned int hit_id_2) const;
	std::pair<double,TMatrixD*> single_jacobian_with_arclength(const genfit::Track& track, const unsigned int hit_id) const;
	std::multimap<double,TMatrixD*> jacobians_with_arclength(const genfit::Track* track) const;
	TVector3 calc_shortest_distance(const TVector3& wire_top, const TVector3& wire_bot, const TVector3& track_pos, const TVector3& track_mom,  TVector3* PCA_on_wire = nullptr, TVector3* PCA_on_track = nullptr) const;
	TRotation calc_rotation_of_vector(const TVector3& v) const;
	TMatrixD rot_to_matrix(const TRotation& rot) const;
	std::vector<TVector3> linear_model_wo_scatter(const genfit::Track& track) const;
	TMatrixD calc_projection_matrix(const TMatrixD& fit_system_base_vectors, const TMatrixD& rotation_global_to_measurement) const;
	void print_model_parameters(const std::vector<TVector3>& model) const;

	bool check_ordered_by_arclen(const genfit::Track& track) const;
	void print_seed_hits(const genfit::Track& track) const;

};

#endif /* CHARMDET_MILLEPEDECALLER_H_ */
