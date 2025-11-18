# MASTER_DMBiosciences_BIOSTAT1

Master Data Management in Biosciences - M1 - Biostatistics I - Université Catholique de Lille

## Course

### STEP 1 – FOUNDATIONS

1.1 – ELEMENTS OF CALCULUS  
1.2 - EPISTEMOLOGY & THEORY OF KNOWLEDGE  

### STEP 2 – STOCHASTIC DYNAMICS & PROBABILITY  

2.1 – MEASURE THEORY  
2.2 – PROBABILITY THEORY  
2.3 – USUAL PROBABILITY DISTRIBUTIONS  
2.4 – ASYMPTOTIC STATISTICS  

### STEP 3 – DATA OBSERVATION

3.1 – DESCRIPTIVE STATISTICS  
3.2 – EXPLORATORY DATA ANALYSIS  

### STEP 4 – INFERENCE & ESTIMATION THEORY

4.1 – PARAMETERS ESTIMATION  
4.2 – EXPERIMENTAL DESIGN & HYPOTHESIS TESTING  
4.3 – TESTING – PARAMETRIC OR NON-PARAMETRIC STATISTICS  
4.4 – BIVARIATE TESTING  
4.5 – PAIRED TESTING  
4.6 – MULTIPLE TESTING CORRECTIONS  

### STEP 5 – LINEAR MODELS EXAMPLES

5.1 – SIMPLE LINEAR REGRESSION  
5.2 – MULTIPLE LINEAR REGRESSION  
5.3 – OTHER REGRESSION MODELS  
5.4 – MODEL SELECTION  

### Ressources 

https://dlsun.github.io/probability/  
https://odr.inrae.fr/intranet/carto/cartowiki/index.php/Statistiques_descriptives_avec_R  
https://chesneau.users.lmno.cnrs.fr/  
https://www.huber.embl.de/msmb/  
https://perso.univ-rennes1.fr/bernard.delyon/regression.pdf  
https://www.sciences.ch/dwnldbl/telecharger.html  
https://www.deeplearningbook.org/contents/linear_algebra.html  
https://www.jstor.org/stable/1390807  
https://cran.r-project.org/doc/contrib/Goulet_introduction_programmation_R.pdf  
https://www.deeplearningbook.org/  
https://pro.univ-lille.fr/fileadmin/user_upload/pages_pros/antoine_ayache/ContinuousRV.pdf  
https://www.datacamp.com/tutorial/probability-mass-function  
https://cran.r-project.org/web/packages/naniar/vignettes/getting-started-w-naniar.html  
https://cran.r-project.org/web/packages/VIM/vignettes/VIM.html  
https://www.r-causal.org/chapters/04-dags#fig-confounder-scatter  
https://delladata.fr/spaghetti-plot-fagot/  
https://www.r-bloggers.com/2025/07/standard-deviation-vs-standard-error-meaning-misuse-and-the-math-behind-the-confusion/  
https://brieflands.com/articles/ijem-71904  
https://academic.oup.com/jrsssc/article/31/2/115/6985178  
https://cran.r-project.org/web/packages/bestNormalize/vignettes/bestNormalize.html  
https://github.com/mboissel/MASTER_OMICS_ANOVA  
https://www.tylervigen.com/spurious-correlations  


### Tips R & co

To start,  

**Explore Packages**,  
•	Browse CRAN (https://cran.r-project.org/), Bioconductor (https://www.bioconductor.org/), and GitHub (https://github.com/) for packages and vignettes.  

**Books & Online References**  
•	The Big Book of R (https://www.bigbookofr.com/) – an extensive compilation of R resources.  
•	R for Data Science (https://r4ds.hadley.nz/) (Hadley Wickham, 2nd ed) – a clear introduction to the tidyverse.  
•	Advanced R (https://adv-r.hadley.nz/) (Hadley Wickham) – deep insights into R’s internals.  

**Websites, Blogs & Forums**  
•	Official CRAN Docs (https://cran.r-project.org/other-docs.html) – guides and tutorials.  
•	frrrenchies (https://frrrenchies.github.io/frrrenchies/) – R tips in French.  
•	R Weekly (https://rweekly.org/) – curated R news.  
•	Jumping Rivers — R Events (https://jumpingrivers.github.io/meetingsR/events.html).  
•	ThinkR (https://thinkr.fr/) – training & articles.  
•	R-bloggers (https://www.r-bloggers.com/) – blog aggregator.  
•	rOpenSci (https://ropensci.org/) – open-source science tools.  
•	CRANberries (https://dirk.eddelbuettel.com/cranberries/) – alerts for new CRAN packages.  
•	French R Forum (CIRAD) (https://forums.cirad.fr/logiciel-R/).  

**Choose your side: Tidyverse (readability) vs data.table (speed)**  
•	[tidyverse inrto](https://juba.github.io/tidyverse/)
•	data.table tutorials:  
  o	[riptutorial.com/data-table](https://riptutorial.com/data-table)
  o	[DataCamp cheat sheet](https://s3.amazonaws.com/assets.datacamp.com/img/blog/data+table+cheat+sheet.pdf)

**Reproducibility & Project Setup**  
•	Ensure reproducibility with renv, Docker, rix, or other tools...  
•	{renv} (dependency management): rstudio.github.io/renv  
•	Docker & Reproducible R:  
  o	[Journey to Reproducibility ](https://m.canouil.dev/journey-reproducibility/#1) 
  o	[Docker_useful repo](https://github.com/mboissel/Docker_useful)  
  o	[Docker feedback PDF ](https://github.com/mboissel/Presentations/blob/master/Docker_feedback_20230718_mboissel.pdf) 
•	{rix} (alternative, rOpenSci): ropensci/rix  
•	Guidance on version control (e.g. ggplot2 versions): brodrigues.co/posts/2025-06-21-ggplot4.html  


**Omics & Biological Data Analysis**  
•	Modern Statistics for Modern Biology (book + labs): huber.embl.de/msmb  
•	Computational Genomics with R: compgenomr.github.io/book  
•	Analysis Script Templates: mboissel/analysis-scripts-templates  
•	DESeq2 experimental design guide: [rstudio-pubs.de…/329027…  ](https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html)
•	mixOmics for multi-omics: mixomics.org  
•	biomartr for annotation: docs.ropensci.org/biomartr

**Advanced & Cross-Language Tools**  
•	consider using VS Code for multi-language projects  
•	reticulate (Python in R): rstudio.github.io/reticulate  
•	useR! Conference on YouTube: [useR! channel  ](https://www.youtube.com/channel/UCv_a9ZGZOH588wUZHZl6T_g/videos)
•	R Conferences: r project.org/conferences  
•	Jenny Bryan on YouTube: [IzRn-QnOhug  ](https://www.youtube.com/watch?v=IzRn-QnOhug&ab_channel=RConsortium)

**Community & Collaboration**  
•	Join R-Ladies or local R User Groups (RUGs) for events and peer learning.  
•	Example: RParis Meetup: meetup.com/rparis  


## Practical Sessions

### Session 1

### Session 2

