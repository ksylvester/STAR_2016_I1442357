// -*- C++ -*-

/// ***************************************************************************************************
/// This is an analysis template created by JW -- changes have been made by KS to attempt an analysis of a specific paper
/// Paper: Jet-Like Correlations with Direct-Photon and Neutral-Pion Triggers at sqrt(sNN)=200 GeV
/// Link to Paper: https://arxiv.org/abs/1604.01117
/// Link to Data: https://drupal.star.bnl.gov/STAR/publications/jet-correlations-direct-photon-and-neutral-pion-triggers-sqrtsnn200gev
/// ***************************************************************************************************


#include "Rivet/HeavyIonAnalysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/ParticleBase.hh"

#define _USE_MATH_DEFINES
#include <math.h>
#include<cmath>
#include<iostream>

namespace Rivet
{

	class STAR_2016_I1442357 : public HeavyIonAnalysis
	{
		public:

			// Constructor
			STAR_2016_I1442357(): HeavyIonAnalysis("STAR_2016_I1442357") {}


			/// Book histograms and initialise projections

			void init()
			{
				HeavyIonAnalysis::init();

				//Select centrality method
				addCentralityMethod(HeavyIonAnalysis::ImpactParameter, 20,"ImpactParameterMethod");

				const FinalState fs;
				addProjection(fs,"FS");

				// Book histograms

				//Figure 1

				_h_CF_1a_gamma = bookHisto1D("d01-x01-y01");
				_h_CF_1a_pi = bookHisto1D("d02-x01-y01");
				_h_CF_1b_gamma = bookHisto1D("d03-x01-y01");
				_h_CF_1b_pi = bookHisto1D("d04-x01-y01");
				_h_CF_1c_gamma = bookHisto1D("d05-x01-y01");
				_h_CF_1c_pi = bookHisto1D("d06-x01-y01");
				_h_CF_1d_gamma = bookHisto1D("d07-x01-y01");
				_h_CF_1d_pi = bookHisto1D("d08-x01-y01");

				//Figure 2

				_h_ZD_2a_au = bookHisto1D("d09-x01-y01");
				_h_ZD_2a_pp = bookHisto1D("d10-x01-y01");
				_h_ZD_2b_au = bookHisto1D("d11-x01-y01");
				_h_ZD_2b_pp = bookHisto1D("d12-x01-y01");

				//Figure 3

				_h_ZD_3_au = bookHisto1D("d13-x01-y01");
				_h_ZD_3_pp = bookHisto1D("d14-x01-y01");

				//Figure 4

				_h_ZD_pp = bookHisto1D("d15-x01-y01");

				//Figure 5

				IAA_photon = bookHisto1D("d16-x01-y01");
				IAA_pi = bookHisto1D("d17-x01-y01");

				//Figure 6

				_h_DZT_photon = bookHisto1D("d18-x01-y01");
				_h_DZT_pi = bookHisto1D("d19-x01-y01");

				//Figure 7

				IAA_trig = bookHisto1D("d20-x01-y01");
				IAA_assoc = bookHisto1D("d21-x01-y01");

				}

			/// Perform the per-event analysis
			void analyze(const Event& event)
			{

				// Get the centrality for each event
				// The first 20 events will give a centrality outside of 0-80 events
				const double c = centrality(event, "ImpactParameterMethod");
				if((c < 0.) || (c > 80.))	vetoEvent;

				const double weight = event.weight();



				// Trigger particle set
				const int pidPI0 = 111;
				const int pidGAMMA = 22;

				// Cuts for trigger particles
				Cuts cutTrigger = Cuts::abseta < 0.9 && ((Cuts::pid == pidPI0 && Cuts::pt > 12.0 * GeV && Cuts::pt < 20.0 * GeV)
				|| (Cuts::pid == pidGAMMA && Cuts::pt > 12.0 * GeV && Cuts::pt < 20.0 * GeV));
				FinalState fs(cutTrigger);
				declare(fs, "partTrigger");

				// Cuts for associated particles
				Cuts cutAssoc = Cuts::abseta < 1.0 && Cuts::pt > 1.2 * GeV && Cuts::pt < 3 * GeV;
				ChargedFinalState cfs(cutAssoc);
				declare(cfs, "partAssoc");


				//Fill Histograms
				double deltaPhi;
				foreach (const Particle& partTrigger, tracksTrigger) {
					// Loop over all associated particles
					foreach (const Particle& partAssoc, tracksAssoc) {
						// Only include associated particles with pT less than the trigger
						if (partAssoc.pt() < partTrigger.pt()) {
							deltaPhi = partAssoc.phi() - partTrigger.phi();

							while (deltaPhi < 0) deltaPhi += 2 * M_PI;

							//Figure 1
							_h_CF_1a_gamma->fill(deltaPhi,1);
							_h_CF_1a_pi->fill(deltaPhi,1);
							_h_CF_1b_gamma->fill(deltaPhi,1);
							_h_CF_1b_pi->fill(deltaPhi,1);
							_h_CF_1c_gamma->fill(deltaPhi,1);
							_h_CF_1c_pi->fill(deltaPhi,1);
							_h_CF_1d_gamma->fill(deltaPhi,1);
							_h_CF_1d_pi->fill(deltaPhi,1);

							//Figure 2

							_h_ZD_2a_au->fill(deltaPhi,1);
							_h_ZD_2a_pp->fill(deltaPhi,1);
							_h_ZD_2b_au->fill(deltaPhi,1);
							_h_ZD_2b_pp->fill(deltaPhi,1);

							//Figure 3

							_h_ZD_3_au->fill(deltaPhi,1);
							_h_ZD_3_pp->fill(deltaPhi,1);

							//Figure 4

							_h_ZD_pp->fill(deltaPhi,1);

							//Figure 5

							IAA_photon->fill(deltaPhi,1);
							IAA_pi->fill(deltaPhi,1);

							//Figure 6

							_h_DZT_photon->fill(deltaPhi,1);
							_h_DZT_pi->fill(deltaPhi,1);

							//Figure 7

							IAA_trig->fill(deltaPhi,1);
							IAA_assoc->fill(deltaPhi,1);
						}
					}
				}
			}	//End Analyze Section

			void finalize()
			{

			}

			/// Name Histograms

			//Figure1

			Histo1DPtr _h_CF_1a_gamma;
			Histo1DPtr _h_CF_1a_pi;
			Histo1DPtr _h_CF_1b_gamma;
			Histo1DPtr _h_CF_1b_pi;
			Histo1DPtr _h_CF_1c_gamma;
			Histo1DPtr _h_CF_1c_pi;
			Histo1DPtr _h_CF_1d_gamma;
			Histo1DPtr _h_CF_1d_pi;

			//Figure 2

			Histo1DPtr _h_ZD_2a_au;
			Histo1DPtr _h_ZD_2a_pp;
			Histo1DPtr _h_ZD_2b_au;
			Histo1DPtr _h_ZD_2b_pp;

			//Figure 3

			Histo1DPtr _h_ZD_3_au;
			Histo1DPtr _h_ZD_3_pp;

			//Figure 4

			Histo1DPtr _h_ZD_pp;

			//Figure 5

			Histo1DPtr IAA_photon;
			Histo1DPtr IAA_pi;

			//Figure 6

			Histo1DPtr _h_DZT_photon;
			Histo1DPtr _h_DZT_pi;

			//Figure 7

			Histo1DPtr IAA_trig;
			Histo1DPtr IAA_assoc;

			//Counters
			int count1=0;
			int count2=0;
			int count3=0;
			int count4=0;

			int Npartt1 =0;
			int Npartp1 =0;
			int Npartt2 =0;
			int Npartp2 =0;
			int Npartt3 =0;
			int Npartp3 =0;
			int Npartt4 =0;
			int Npartp4 =0;

	};

	// The hook for the plugin system
	DECLARE_RIVET_PLUGIN(STAR_2016_I1442357);

}
