//! Case-level information for use in models.

use crate::Time;

/// Information generated about a case.
///
/// This type is only relevant when implementing disease models.
///
/// Times should be provided relative to the exposure time of the case.
#[derive(Debug)]
pub struct CaseHistory {
    /// Defines how the infectiousness of the case changes over time.
    ///
    /// It is a indexed by time steps relative to the onset of infectiousness, so
    /// `infectiousness[0]` gives the infectiousness at the time of onset (end of incubation
    /// period) and `infectiousness[1]` is the infectiousness one time step later etc.
    ///
    /// The total area under the described curve (i.e. the sum of elements) is equal to the
    /// expected number of new individuals that this case will infect.
    pub infectivity: Vec<f64>,

    /// Time after exposure when symptom onset occurs.
    pub symptom_onset: Option<Time>,

    /// Time after exposure when the case is reported.
    pub reported: Option<Time>,
}

#[derive(Debug)]
pub(crate) enum Case {
    Latent(Vec<f64>),
    Active(Vec<f64>),
    Recovered,
}

/// Important events in the disease history of a case.
///
/// These are times relative to the start of the outbreak.
#[derive(Debug)]
pub struct History {
    /// Time when a case was initially infected.
    pub infected: Time,

    /// Time when the case became infectious.
    pub infectious_onset: Time,

    /// Time when the case was most infectious.
    pub infectious_peak: Time,

    /// Time when the case recovered.
    pub recovered: Time,

    /// Time when the case would have been reported.
    pub reported: Option<Time>,

    /// Time when symptoms began, if at all.
    pub symptom_onset: Option<Time>,
}

impl CaseHistory {
    pub(crate) fn into_case_history(self) -> (Case, History) {
        let milestones = milestones(&self.infectivity);

        (
            Case::Latent(self.infectivity.into_iter().rev().collect()),
            History {
                infected: 0,
                infectious_onset: milestones[0],
                infectious_peak: milestones[1],
                recovered: milestones[2],
                reported: self.reported,
                symptom_onset: self.symptom_onset,
            },
        )
    }
}

fn milestones(infectivity: &[f64]) -> [Time; 3] {
    let mut ever_infectious = false;
    let mut onset = 0;
    let (mut peak, mut peak_val) = (0, 0.0);
    let mut recov = 0;

    for (t, &inf) in infectivity.iter().enumerate() {
        let t = Time::try_from(t).unwrap();

        if inf > peak_val {
            (peak, peak_val) = (t, inf);
        }

        if !ever_infectious && inf > 0.0 {
            ever_infectious = true;
            onset = t;
        }

        if ever_infectious && inf <= 0.0 {
            recov = t;
            break;
        }

        recov = t + 1;
    }

    [onset, peak, recov]
}

impl Case {
    pub(crate) fn is_recovered(&self) -> bool {
        match self {
            Case::Latent(_) | Case::Active(_) => false,
            Case::Recovered => true,
        }
    }

    pub(crate) fn step(&mut self) -> f64 {
        match self {
            Case::Latent(inf) => match inf.pop() {
                None => {
                    *self = Case::Recovered;
                    0.0
                }
                Some(i) if i > 0.0 => {
                    *self = Case::Active(inf.clone());
                    i
                }
                Some(i) => i,
            },
            Case::Active(inf) => match inf.pop() {
                None => {
                    *self = Case::Recovered;
                    0.0
                }
                Some(i) if i <= 0.0 => {
                    *self = Case::Recovered;
                    0.0
                }
                Some(i) => i,
            },
            Case::Recovered => 0.0,
        }
    }
}

impl History {
    pub(crate) fn time_shift_forward(&mut self, offset: Time) {
        self.infected += offset;
        self.infectious_onset += offset;
        self.infectious_peak += offset;
        self.recovered += offset;
        for time in &mut self.reported {
            *time += offset;
        }
        for time in &mut self.symptom_onset {
            *time += offset;
        }
    }

    pub(crate) fn time_shift_back(&mut self, offset: Time) {
        self.infected -= offset;
        self.infectious_onset -= offset;
        self.infectious_peak -= offset;
        self.recovered -= offset;
        for time in &mut self.reported {
            *time -= offset;
        }
        for time in &mut self.symptom_onset {
            *time -= offset;
        }
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = Time> {
        [
            Some(self.infected),
            Some(self.infectious_onset),
            Some(self.infectious_peak),
            Some(self.recovered),
            self.symptom_onset,
            self.reported,
        ]
        .into_iter()
        .flatten()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_milestones() {
        let inf = &[0.0, 0.0, 0.1, 0.2, 0.6, 0.1, 0.0];
        assert_eq!(milestones(inf), [2, 4, 6]);

        let inf = &[0.0, 0.0, 0.1, 0.2, 0.6, 0.1];
        assert_eq!(milestones(inf), [2, 4, 6]);

        let inf = &[0.0, 0.0, 0.1, 0.2, 0.6, 0.1, 0.0, 1.0];
        assert_eq!(milestones(inf), [2, 4, 6]);

        let inf = &[0.1, 0.2, 0.6, 0.1, 0.0, 1.0];
        assert_eq!(milestones(inf), [0, 2, 4]);
    }
}
