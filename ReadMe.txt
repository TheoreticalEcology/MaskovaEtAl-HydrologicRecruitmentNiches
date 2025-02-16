Original germination data are in table "data.csv". Each row in this datatable represents one pot in germination experiment. Column A represent from which waterbody the sediment was collected, column B represent the hydrological treatment during the germination experiment, columns C - DD represents the number of seedlings per respective plant species which germinated during the experiment. 

Additional species-specific data collected in literature are in table "traits.csv".

Scripts "species_specific_permanent.R" and "species_specific_temporal.R" - fitted Kruskal-Wallis test for individual species and permanent and temporal hydrological treatments respectively. Post-hoc dunn's test was used to compare each hydrological treatment pairwise. Results were plotted as heatmap figures.


Scripts "models_permanent.R" and "models_temporal.R" - fitted three separate models, each with one seed trait (dormancy, fluctuating temperature requirements and oxygen requirements for germination) and its interaction with hydfological treatments. Model residuals, including a test of phylogenetic signal in the residuals based on Moran's I test, using the DHARMa package.
Results were plotted as heatmap figures.

