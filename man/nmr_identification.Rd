\name{nmr_identification}
\alias{nmr_identification}
\title{NMR metabolite identification}
\description{
This function performs metabolite identification on a dataset of nmr peaks.
}
\usage{
nmr_identification(dataset, ppm.tol, frequency_scores, solvent_scores, organism_scores,
                   method='Match_uniq', per.sample=FALSE, tresh_zero=0, alpha=10e-4)
}
\arguments{
\item{dataset}{List representing the dataset from an nmr peaks metabolomics experiment.}

\item{ppm.tol}{ppm tolerance when matching reference peaks to the dataset peaks.}

\item{per.sample}{Logical indicating wether identification should be done for each sample (for each sample, only the peaks present in that sample will be used in the identification) or for all samples together (all peaks will be used together in the identification).}

\item{tresh_zero}{Intensity value above which peaks are considered present. Only necessary if per.sample=TRUE}

\item{method}{Identification method to use. There are three: 'Match_uniq' (default), 'Hyper' or 'Hyper_uniq'. For more information, see details below.}

\item{alpha}{A metabolite with a corrected p.value above alpha will not be considered identified. Only necessary if method is either 'Hyper' or 'Hyper_uniq'.}

\item{frequency_scores}{Each frequency in the library should be given a score from 0 to 1, according to the frequency in which the samples were acquired. A score of 1 should be given to the frequency under which the samples were acquired. Argument should be in the form of a list, like: list('400'=0, '500'=1, '600'=0, '700'=0). For more information, see details below.}

\item{solvent_scores}{Each solvent in the library should be given a score from 0 to 1, according to the solvent with which the samples were acquired. A score of 1 should be given to the solvent with which the samples were acquired. Argument should be in the form of a list, like: list(CD3OD=0, D2O=1, Water=1, CDCl3=0, 'Acetone-d6'=0, 'Acetone'=0, "DMSO-d6"=0, "100\%_DMSO"=0, "5\%_DMSO"=0, C=0, C6D6=0, CD3CN=0, C2D2Cl4=0, CD2Cl2=0, CDC3OD=0, Ethanol=0). For more information, see details below.}

\item{organism_scores}{This gives the opportunity to score each reference according to the presence or not of the respective metabolite in the organism(s) or group(s) of organisms under study, so that metabolites not present in an organism/group of interest will have a lower score and thus be less likely to be identified. The presence of a metabolite in an organism/group is evaluated using the information in the KEGG database. Thus, all organism/group codes should be present in KEGG. Argument should be in the form of a list, like: list('hsa'=1, 'other'=0, 'not_in_kegg'=0). For more information, see details below.}
}
\value{
If per.sample=FALSE, it gives a list with two items: results_table and more_results. If per.sample=TRUE, it gives a list with as many items as the number of samples. Each sample contains a list with the two items results_table and more_results.

results_table is a data.frame with the information on the metabolites matched. Each row corresponds to a spectrum from the library that matched the dataset:
\describe{
   \item{SPCMNM}{ID of the metabolite in our library.}
   \item{Name}{Name of the metabolite.}
   \item{SPCMNS}{ID of the respective spectrum in our library.}
   \item{Final_Score}{Final score.}
   \item{match_score}{Matching score.}
   \item{hypergeometric_score}{Hypergeometric score, if method='Hyper' or 'Hyper_uniq'.}
   \item{ratio}{matched ratio score, if method='Match_uniq'}
   \item{uniqueness_score}{uniqueness_score, if method = 'Hyper_uniq' or 'Match_uniq'}
   \item{score_frequency}{Score given to the frequency of that spectrum.}
   \item{score_solvents}{Score given to the solvent of that spectrum.}
   \item{score_organisms}{The organism score, according to the metabolite's presence in one of the organisms/groups given by the user in organism_scores argument.}
   \item{n.peaks.matched}{Number of peaks from the metabolite's spectrum that matched the sample.}
   \item{detailed_results_id}{ID to access the more detailed results in the item more_results.}
}
   
more_results is a list whose items are identified by an ID that is specified in the detailed_results_ID column of the results_table. Each item is a list with the following information:
\describe{
    \item{matched_peaks_ref}{Vector with the peaks from the reference spectrum that matched the sample. Reference peak in the ith position matched the sample peak in the ith position of the vector matched_peaks_samp.}
    \item{matched_peaks_samp}{Vector with the peaks from the sample that matched the reference spectrum. Sample peak in the ith position matched the reference peak in the ith position of the vector matched_peaks_ref.}
    \item{reference_peaks}{Vector with all the peaks in the reference spectrum.}
}
}
\details{
There are three methods implemented to perform metabolite identification. The default one is Matched Ratio with uniqueness score ('Match_uniq'):

- Hypergeometric Test: it calculates the probability of a group of k peaks matching to a certain reference spectrum not being caused by chance. A p-value over alpha denotes that the metabolite corresponding to the reference spectrum in question is not present in the sample. After all reference metabolites are matched to the samples, the p-values are adjusted for multiple testing using the False Discovery Rate (FDR) method. For those with p.value under alpha, the score is transformed into a scale of 0 to 1, by applying the following: 1-(p.value/alpha).

- Hypergeometric Test with Uniqueness score ('Hyper_uniq'): the final score of a reference spectrum is obtained by calculating the average between the hypergeometric test score and the uniqueness score, if the hypergeometric test score is not null. The uniqueness score of a reference is the average of the uniqueness rate of all peaks in that reference. The uniqueness rate of a peak is calculated by dividing 1 by the number of reference spectra that peak is in, based on the reference library used.

- Matched ratio scores with uniqueness score ('Match_uniq'): the final score of a reference spectrum is obtained by calculating the average between the matched ratio score and the uniqueness score, if the matched ratio score is not null. The matched ratio score gives the ratio of peaks from a reference that matched the sample. It is calculated by dividing the number of different peaks matched between the reference and the sample by the total number of different peaks in such reference.

After scoring a match between a reference and a sample, the reference is further scored regarding the conditions under which it was acquired. To do so, each frequency and solvent represented in the library must be given a score by the user. The user also has the possibility to score each reference according to the presence or not of the respective metabolite in the organism(s) or group(s) of organisms under study, so that metabolites not present in an organism/group of interest will have a lower score and thus be less likely to be identified. The presence of a metabolite in an organism/group is evaluated using the information in the KEGG database.
}
\examples{
\donttest{
 propolis_mv=missingvalues_imputation(propolis)
 freq_scores = list('400'=0, '500'=0, '600'=1, '700'=0)
 solv_scores = list(CD3OD=0, D2O=.8, Water=.8, CDCl3=0, 'Acetone-d6'=0, 'Acetone'=0,"DMSO-d6"=0,
                    '100\%_DMSO'=0,
                    '5\%_DMSO'=0, C=0, C6D6=0, CD3CN=0, C2D2Cl4=0, CD2Cl2=0,
                    CDC3OD=1, Ethanol=0)
 org_scores = list('Eudicots'=1, 'Monocots'=1, 'ame'=.9, 'other'=0, 'not_in_kegg'=0)
 id_res=nmr_identification(propolis_mv, ppm.tol=0.03, freq_scores, solv_scores, org_scores)
}
}
