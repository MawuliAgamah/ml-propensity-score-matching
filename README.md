# Estimating intervention effect’s using propensity score matching and artificial intelligence
 
<p> The files are structured as follows: </p>

<p> Jupyter notebooks -> lalonde.ipynb - Exploratory data analysis and ML for propensity score estimation (Logistic regression,classification tree,Random forest,extreme gradident boosted tree ,artificial neurual netowork) </p>

<p>The seperate artificial neural network notebooks implement ANN's using pytorch alone, however there are issues with model's losses,the model used in our thesis is implemented with skorch (in the main lalonde notebook), this allows allows us to use sklearns grid search. </p>

<p> Jupyter notebooks -> syntheticstudies.ipynb - simulated experimental and nonexperiemntal datasets for testing models given ideal conditions (large sample sizes, balanced covariates) </p> 
<p> R-script -> Matched samplea nlysis.R - Matching , stratification and average treatment effect estimates  </p>
