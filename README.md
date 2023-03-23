---
editor_options: 
  markdown: 
    wrap: 72
---

Description of code and data files affiliated with **Social environment
influences the temporal dynamics of sneak-spawning in a fish with
alternative reproductive tactics.**

[**Authors**]{.ul}:

Matthew C. Kustra\*, Kelly A. Stiver, Susan Marsh-Rollo, Jennifer K.
Hellmann, and Suzanne H. Alonzo

\*For questions
contact:[mkustra\@ucsc.edu](mailto:mkustra@ucsc.edu){.email}

[**Abstract**]{.ul}:

Several predictions of sperm competition theory are not well supported
empirically. One potential reason is that most current theory and
empirical research ignore how the social environment influence the
temporal dynamics of mating. We propose that understanding these
dynamics is key to understanding sexual selection and improving the
predictive power of theory. To demonstrate the importance of these
dynamics, we quantify how males' social role, interactions among males,
and current social environment influence the timing of mating in
Symphodus ocellatus, a species with three alternative male reproductive
tactics. Nesting males spawn synchronously with females; sneakers and
satellites sneak-spawn with some time-delay. Satellites also cooperate
with nesting males. We found that satellites have shorter sneak-spawning
delays than sneakers, a benefit of their cooperation with nesting males.
Sneak-spawning delays decreased with increasing nest activity for
sneakers but not satellites, suggesting sneakers may benefit from
increased sperm competition intensity. Current sperm competition models
ignore this potential benefit which may be why the prediction that males
should decrease investment when sperm competition involves more than two
males is not well supported. Our study provides insight into mechanisms
that drive variation in the timing of spawning, which could explain
mismatches between theoretical and empirical results.

[**Code**]{.ul}: To properly run code make sure "**Data.csv**" is
placed in a folder named "**Data**" and all model results (files ending in .rds) are in a folder named **ModelResults**. Data files are on [dryad](https://doi.org/10.7291/D17698).

-   **"Analysis.R":** R script file that runs all statistical analyses
    and makes all the figures from the paper (excluding figure 1). Uses
    "**Data.csv"** in the "**Data**" folder. Code is available on Zenodo
    and
    [github](https://github.com/Kustra-Matt/SocialEnvironmentAffectsSneakTimes)

-   **"sessioninfo.txt":** text file with the session information for
    when these analyses were performed(version of R, package versisons
    used, computer hardware, etc.).

[**Data:**]{.ul} Folder with "**Data.csv**"

-   "**Data.csv**": CSV file that contains all the data used for this
    paper. Below are a description of each column:

    -   ***Nest_ID*****:** Name of the nest.

    -   ***File.Name*****:** Name of the video file for the observation.

    -   ***Lag.Time*****:** Time delay (seconds) from when nesting male
        spawns and sneaker/satellite male reaches the nest.

    -   ***Male*****:** Male type. SN = sneaker; SAT = Satellite male.

    -   ***SpawnNumber*****:** Number of spawning event in order of
        occurrence for a given video file.

    -   ***SN_Sneak_S*****:** Number of sneaker that sneaked in a given
        spawning event. Will be the same number for all observations
        with the same file.name/nest_ID and spawn number.

    -   ***Sat_Sneak_S*****:** Number of satellite sneaks in a given
        spawning event (0 or 1). Will be the same number for all
        observations with the same file.name/nest_ID and spawn number.

    -   ***TotalSneaks_S*****:** Total number of males that sneaked in a
        given spawning event (**SN_Sneak_S** + **Sat_Sneak_S**). Will be
        the same number for all observations with the same
        file.name/nest_ID and spawn number.

    -   ***AloneS*****:** Category of spawning event: "Alone" = either
        one sneaker or one satellite male, "Paired" = one sneaker and
        one satellite, and "Many" = one satellite and multiple sneakers.
        Will be the same number for all observations with the same
        file.name/nest_ID and spawn number.

    -   ***SN_Sneak_N*****:** Total number of sneaks performed by a
        sneaker at a given nest. Will be the same number for all
        observations from the same file.name/nest_ID.

    -   ***SAT_Sneak_N*****:** Total number of sneaks performed by a
        satellite male at a given nest. Will be the same number for all
        observations from the same file.name/nest_ID.

    -   ***Sneak_N*****:** Total number of sneaks at a given nest
        (**SN_Sneak_N** + **SAT_Sneak_N).** Will be the same number for
        all observations from the same file.name/nest_ID.

    -   ***TOTAL.FEMALE.VISITS*****:** Total number of females that
        visited the nest during the ten minute observation*.* Will be
        the same number for all observations from the same
        file.name/nest_ID.

    -   ***TOTAL.FEMALES.SPAWNING*****:** Total number of females that
        spawned at the nest during the ten minute observation. Will be
        the same number for all observations from the same
        file.name/nest_ID.

    -   ***TOTAL.SPAWNS*****:** Total number of spawns at the nest
        during the ten minute observation. Will be the same number for
        all observations from the same file.name/nest_ID.

    -   ***Nest.male.to.SAT_AG*****:** Number of aggressive behaviors
        from nesting male to satellite. Will be the same number for all
        observations from the same file.name/nest_ID.

    -   ***SAT.to.SN_AG*****:** Number of aggressive behaviors from
        satellite male to sneakers.

    -   ***Sat.to.NM_SUB*****:** Number of submissive behaviors from
        satellite male to nesting male. Will be the same number for all
        observations from the same file.name/nest_ID.

    -   ***Nest.male.to.SN_AG*****:** Number of aggressive behaviors
        from nesting male to sneakers. Will be the same number for all
        observations from the same file.name/nest_ID.

    -   ***Average.Sn*****:** Average number of sneakers during the ten
        minute observation. Will be the same number for all observations
        from the same file.name/nest_ID.

[**ModelResults:**]{.ul} Folder with rds files of statistical model
outputs. These are not raw data and are not needed to replicate any of
the analyses, just the results of the analyses in the code. They are
provided to make replication of graphs faster because these statistical
models can take a long time to run.

-   "**PCA_Randomization_Results.rds**": R data file that contains the
    PCA randomization result output. Used to make SI Figures 2-4. Since
    this analysis can take a while one can just load up the results for
    faster visualization.

-   "**Model_1.rds**": R data file that contains the Bayesian model
    result output for the first analysis of sneaker and satellite time
    delays (results reported in Table 1 and Figure 2). Since this
    analysis can take a while one can just load up the results for
    faster visualization.

-   "**Model_SN_Only.rds**":R data file that contains the Bayesian model
    result output for the analysis of sneaker only spawns (results
    reported for Table and satellite time delays (results reported in SI
    Table 4 and SI Figure 5). Since this analysis can take a while one
    can just load up the results for faster visualization.

-   "**Model_1\_Regularized_Prior.rds**": R data file that contains the
    Bayesian model result output for the first analysis sensitivity test
    of using a regularizing prior (results reported in SI Table1). Since
    this analysis can take a while one can just load up the results for
    faster visualization."

-   "**Model_2.rds**": R data file that contains the Bayesian model
    result output for the second analysis of sneaker and satellite time
    delays using the nest-level effects (results reported in Table 2 and
    Figure 3). Since this analysis can take a while one can just load up
    the results for faster visualization.

-   "**Model_2\_Regularized_Prior.rds**": R data file that contains the
    Bayesian model result output for the second analysis sensitivity
    test of using a regularizing prior (results reported in SI Table2).
    Since this analysis can take a while one can just load up the
    results for faster visualization.

-   "**Model_3.rds**": R data file that contains the Bayesian model
    result output for the third analysis that included nest and spawn
    level effects (results reported in SI Table3). Since this analysis
    can take a while one can just load up the results for faster
    visualization.
