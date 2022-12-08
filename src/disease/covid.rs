//! Disease model tuned for Covid-19.

use super::DiseaseModel;
use crate::case::CaseHistory;
use crate::Time;
use rand::seq::index::sample_weighted;
use rand::Rng;
use rand_distr::Distribution;

/// Covid-19 disease model.
///
/// From exposure there is a 2 day latent period before the onset of infectiousness. This is
/// followed by a 3 day exponential increase to peak infectiousness, then a linear decline over the
/// following 12 days.
///
/// When cases do not have symptoms, their infectivity is multiplied by a factor of 0.3. For cases
/// that do develop symptoms, the onset will occur within the first 3 days and last until the end
/// of infectiousness.
///
/// Based on the model used in
/// Chang, S.L., Harding, N., Zachreson, C. et al.
/// Modelling transmission and control of the COVID-19 pandemic in Australia.
/// _Nat Commun_ **11**, 5710 (2020).
/// doi: [10.1038/s41467-020-19393-6](https://doi.org/10.1038/s41467-020-19393-6)
#[derive(Debug)]
pub struct Covid<DRep, DR> {
    /// A `Distribution<Time>` of times between symptom onset and case notification.
    pub reporting_time: DRep,

    /// A `Distribution<f64>` of individual reproduction numbers.
    ///
    /// Note that this applies to the average number of cases infected in the baseline, but any
    /// individuals who are asymptomatic or pre-symptomatic will have a reduced chance of
    /// transmission for the relevant period. This means that the actual number of cases will tend
    /// to be fewer than implied by this distribution alone.
    pub reproduction_number: DR,
}

/// Baseline infectiousness curve normalised to unity.
#[rustfmt::skip]
const INFECTIOUSNESS: &[f64] = &[
    0.0, 0.0, // latent
    0.04, 0.08, 0.16, // exponential growth
    0.144, 0.128, 0.112, 0.096, 0.08, 0.064, 0.048, 0.032, 0.016, // linear drop
];

/// Fraction of baseline infectiousness when asymptomatic/pre-symptomatic.
const ASYMP_INFECT: f64 = 0.3;

/// Proportion of cases that will ever become symptomatic.
const FRAC_SYMPTOMATIC: f64 = 0.667;

/// Relative proportion of symptomatic cases that will have symptom onset on each day.
///
/// i.e. 30% develop symptoms on the day of infection, 50% on the following day.
const SYMPTOMS: [u32; 3] = [3, 5, 2];

impl<DRep, DR> DiseaseModel for Covid<DRep, DR>
where
    DRep: Distribution<Time>,
    DR: Distribution<f64>,
{
    type State = ();

    fn generate_case<R: Rng>(&self, _state: &mut Self::State, mut rng: R) -> CaseHistory {
        let will_have_symptoms = rng.gen_bool(FRAC_SYMPTOMATIC);
        let symptom_onset: Option<Time> = will_have_symptoms.then(|| {
            sample_weighted(&mut rng, SYMPTOMS.len(), |x| SYMPTOMS[x], 1)
                .unwrap()
                .index(0) as Time
        });
        let reported = symptom_onset.map(|x| x + self.reporting_time.sample(&mut rng));

        let mut infectivity = INFECTIOUSNESS.to_vec();

        // reduce infectivity for pre-symptomatic period or is asymptomatic
        let end_pos = symptom_onset
            .map(|x| x as usize)
            .unwrap_or(INFECTIOUSNESS.len());
        for inf in &mut infectivity[0..end_pos] {
            *inf *= ASYMP_INFECT;
        }

        // scale by reproduction number
        let r = self.reproduction_number.sample(&mut rng);
        for inf in &mut infectivity {
            *inf *= r;
        }

        CaseHistory {
            infectivity,
            symptom_onset,
            reported,
        }
    }

    // fn generate_singleton<R: Rng>(&self, mut rng: R) -> CaseHistory {
    // }
}
