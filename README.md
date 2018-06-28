<h1><p align="center">
Manual: uCAREChemSuiteCLI
</p></h1>

-   [Introduction](#introduction)
-   [Requirements](#requirements)
-   [Imports](#imports)
-   [Installation](#installation)
-   [Example](#example)
-   [Usage](#usage)
-   [Help](#help)
-   [Bugs](#bugs)
-   [License](#license)

Introduction
============
<p>This package consist of four functions viz. "drug.class.deterministic", "drug.class.stochastic", "drug.resistome.deterministic" and "drug.resistome.stochastic" respectively to predict resistome of <i>Escherichia coli</i> for drug chemical structure.</p>
<p>
It predicts the class/family of unknown candidate drug molecule using two algorithms viz. deterministic model ("drug.class.deterministic") and stochastic model ("drug.class.stochastic") [Unpublished data]. Furthermore once the drug class is predicted, the resistome of the predicted drug class can be fetched out from database using "drug.resistome.deterministic" and "drug.resistome.resistome" functions.</p>

Requirements
============
-   R (tested in version 3.4.3)

Imports
============
-   ChemmineR (tested in version 2.30.2)
-   stats (tested in version 3.4.3)
-   utils (tested in version 3.4.3)
-   usethis (tested in version 1.3.0)

Installation
============
The package can be installed from [CRAN](https://cran.r-project.org/package=uCAREChemSuiteCLI) (recommended). 

```R
## CMD Installation
 install.packages("uCAREChemSuiteCLI")
```

Example
============
We have provided examplary file under ```uCAREChemSuiteCLI/inst/extdata/```

Usage
=====
Run the algorithm using the provided examples with the following command:

```R
# Loading libraries
 library("uCAREChemSuiteCLI")
```

```R
# Run drug.class.deterministic
example.class.deterministic<- system.file('extdata/example.sdf', package="uCAREChemSuiteCLI")
drug.class.deterministic(example.class.deterministic)
```

```R
# Run drug.class.stochastic
example.class.stochastic<- system.file('extdata/example.sdf', package="uCAREChemSuiteCLI")
drug.class.stochastic(example.class.stochastic,"3","0.25")
```

```R
# Run drug.resistome.deterministic
example.resistome.deterministic<- system.file('extdata/example.sdf', package="uCAREChemSuiteCLI")
drug.resistome.deterministic(example.resistome.deterministic)
```

```R
# Run drug.resistome.stochastic
example.resistome.stochastic<- system.file('extdata/example.sdf', package="uCAREChemSuiteCLI")
drug.resistome.stochastic(example.resistome.stochastic, "3", "0.25")
```

Help
============
All functions are documented. You can find additional information using the help function of R. 
<br> eg. `??drug.class.deterministic`

Bugs
===========
Please report any issues or bugs you find while installing/running **uCAREChemSuiteCLI** to:
-   Saurav Bhaskar Saha [<saurav.saha@shiats.edu.in>]

License
============
uCAREChemSuiteCLI is under MIT License.
