How to load and use data
------------------------

This package contains data in the SummarizedExperiment formant and saved
as .Rds files. This data structure contains all assays and metadata.
Files can be opened with R, which is a free software. First install R on
your system (we recomend through
<https://www.rstudio.com/products/rstudio/>)

next you will have to install the SummarizedExperiment library in order
to efficiently access the data. for this type

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("SummarizedExperiment")

Once installed the package can be loaded using

    library(SummarizedExperiment)

Now, individual files can be loaded using:

    data <- readRDS("seq/TimeSeriesFST.dHC.SE.rds")

You can access the count data by using:

    countdat <- assays(data)$counts

lets inspect the data (top-left 4x4 sub-matrix)

    countdat[1:4,1:4]

    ##               TimeSeriesFST_id_1.dHC.Homecage.0min
    ## 0610009O20Rik                                  661
    ## 1700028P14Rik                                   17
    ## 1810055G02Rik                                  371
    ## 2310039H08Rik                                   32
    ##               TimeSeriesFST_id_12.dHC.Homecage.0min
    ## 0610009O20Rik                                   349
    ## 1700028P14Rik                                    16
    ## 1810055G02Rik                                   199
    ## 2310039H08Rik                                     3
    ##               TimeSeriesFST_id_13.dHC.Homecage.0min
    ## 0610009O20Rik                                   632
    ## 1700028P14Rik                                    14
    ## 1810055G02Rik                                   296
    ## 2310039H08Rik                                    24
    ##               TimeSeriesFST_id_22.dHC.Homecage.0min
    ## 0610009O20Rik                                   746
    ## 1700028P14Rik                                    21
    ## 1810055G02Rik                                   355
    ## 2310039H08Rik                                    45

All metadata can be accessed by using:

    colData(data)

    ## DataFrame with 43 rows and 11 columns
    ##                                           SampleNumber        AnimalNumber
    ##                                               <factor>            <factor>
    ## TimeSeriesFST_id_1.dHC.Homecage.0min  TimeSeriesFST_2  TimeSeriesFST_id_1 
    ## TimeSeriesFST_id_12.dHC.Homecage.0min TimeSeriesFST_24 TimeSeriesFST_id_12
    ## TimeSeriesFST_id_13.dHC.Homecage.0min TimeSeriesFST_26 TimeSeriesFST_id_13
    ## TimeSeriesFST_id_22.dHC.Homecage.0min TimeSeriesFST_44 TimeSeriesFST_id_22
    ## TimeSeriesFST_id_25.dHC.Homecage.0min TimeSeriesFST_50 TimeSeriesFST_id_25
    ## ...                                                ...                 ...
    ## TimeSeriesFST_id_27.dHC.Swim.4h       TimeSeriesFST_54 TimeSeriesFST_id_27
    ## TimeSeriesFST_id_3.dHC.Swim.4h        TimeSeriesFST_6  TimeSeriesFST_id_3 
    ## TimeSeriesFST_id_31.dHC.Swim.4h       TimeSeriesFST_62 TimeSeriesFST_id_31
    ## TimeSeriesFST_id_40.dHC.Swim.4h       TimeSeriesFST_80 TimeSeriesFST_id_40
    ## TimeSeriesFST_id_7.dHC.Swim.4h        TimeSeriesFST_14 TimeSeriesFST_id_7 
    ##                                         Region Condition TimePoint
    ##                                       <factor>  <factor>  <factor>
    ## TimeSeriesFST_id_1.dHC.Homecage.0min       dHC  Homecage      0min
    ## TimeSeriesFST_id_12.dHC.Homecage.0min      dHC  Homecage      0min
    ## TimeSeriesFST_id_13.dHC.Homecage.0min      dHC  Homecage      0min
    ## TimeSeriesFST_id_22.dHC.Homecage.0min      dHC  Homecage      0min
    ## TimeSeriesFST_id_25.dHC.Homecage.0min      dHC  Homecage      0min
    ## ...                                        ...       ...       ...
    ## TimeSeriesFST_id_27.dHC.Swim.4h            dHC      Swim        4h
    ## TimeSeriesFST_id_3.dHC.Swim.4h             dHC      Swim        4h
    ## TimeSeriesFST_id_31.dHC.Swim.4h            dHC      Swim        4h
    ## TimeSeriesFST_id_40.dHC.Swim.4h            dHC      Swim        4h
    ## TimeSeriesFST_id_7.dHC.Swim.4h             dHC      Swim        4h
    ##                                                      Block CageNumber      Sex
    ##                                                   <factor>  <integer> <factor>
    ## TimeSeriesFST_id_1.dHC.Homecage.0min  TimeSeriesFST_Block1     227127     Male
    ## TimeSeriesFST_id_12.dHC.Homecage.0min TimeSeriesFST_Block2     227129     Male
    ## TimeSeriesFST_id_13.dHC.Homecage.0min TimeSeriesFST_Block3     227130     Male
    ## TimeSeriesFST_id_22.dHC.Homecage.0min TimeSeriesFST_Block4     227132     Male
    ## TimeSeriesFST_id_25.dHC.Homecage.0min TimeSeriesFST_Block5     227133     Male
    ## ...                                                    ...        ...      ...
    ## TimeSeriesFST_id_27.dHC.Swim.4h       TimeSeriesFST_Block5     227133     Male
    ## TimeSeriesFST_id_3.dHC.Swim.4h        TimeSeriesFST_Block1     227127     Male
    ## TimeSeriesFST_id_31.dHC.Swim.4h       TimeSeriesFST_Block6     227134     Male
    ## TimeSeriesFST_id_40.dHC.Swim.4h       TimeSeriesFST_Block7     227136     Male
    ## TimeSeriesFST_id_7.dHC.Swim.4h        TimeSeriesFST_Block2     227128     Male
    ##                                          Experiment SampleProcessing Condition2
    ##                                            <factor>         <factor>   <factor>
    ## TimeSeriesFST_id_1.dHC.Homecage.0min  TimeSeriesFST         WholeRNA   Homecage
    ## TimeSeriesFST_id_12.dHC.Homecage.0min TimeSeriesFST         WholeRNA   Homecage
    ## TimeSeriesFST_id_13.dHC.Homecage.0min TimeSeriesFST         WholeRNA   Homecage
    ## TimeSeriesFST_id_22.dHC.Homecage.0min TimeSeriesFST         WholeRNA   Homecage
    ## TimeSeriesFST_id_25.dHC.Homecage.0min TimeSeriesFST         WholeRNA   Homecage
    ## ...                                             ...              ...        ...
    ## TimeSeriesFST_id_27.dHC.Swim.4h       TimeSeriesFST         WholeRNA    Swim 4h
    ## TimeSeriesFST_id_3.dHC.Swim.4h        TimeSeriesFST         WholeRNA    Swim 4h
    ## TimeSeriesFST_id_31.dHC.Swim.4h       TimeSeriesFST         WholeRNA    Swim 4h
    ## TimeSeriesFST_id_40.dHC.Swim.4h       TimeSeriesFST         WholeRNA    Swim 4h
    ## TimeSeriesFST_id_7.dHC.Swim.4h        TimeSeriesFST         WholeRNA    Swim 4h

To see all the different types of processed data you can use

    names(assays(data))

    ## [1] "counts" "logcpm" "tpm"

lets access logcpm data instead of count data (top-left 4x4 sub-matrix)

    assays(data)$logcpm[1:4,1:4]

    ##               TimeSeriesFST_id_1.dHC.Homecage.0min
    ## 0610009O20Rik                             3.691646
    ## 1700028P14Rik                             0.696081
    ## 1810055G02Rik                             3.133395
    ## 2310039H08Rik                             1.062437
    ##               TimeSeriesFST_id_12.dHC.Homecage.0min
    ## 0610009O20Rik                             3.6575792
    ## 1700028P14Rik                             1.0048348
    ## 1810055G02Rik                             3.1150689
    ## 2310039H08Rik                             0.2811467
    ##               TimeSeriesFST_id_13.dHC.Homecage.0min
    ## 0610009O20Rik                             3.7534576
    ## 1700028P14Rik                             0.6539037
    ## 1810055G02Rik                             3.0211836
    ## 2310039H08Rik                             0.9486970
    ##               TimeSeriesFST_id_22.dHC.Homecage.0min
    ## 0610009O20Rik                             3.8443178
    ## 1700028P14Rik                             0.8273351
    ## 1810055G02Rik                             3.1250077
    ## 2310039H08Rik                             1.3239717

Certain datasets have further data associated with it. lets have a look
at the phosphodata for example

    data_phos <- readRDS("phos/PhosphoData_RMAmbig.SE.rds")

lets see the variance stabilized data (top-left 4x4 sub-matrix)

    assays(data_phos)$vsn[1:4,1:4]

    ##                                                                                  Exp2_dHC_SW15_1
    ## _AAAAAAAAAAAAGDSDSWDADTFSMEDPVR_                                                        17.05760
    ## _AAAAALSGS[Phospho_STY_NL_98 (STY)]PPQTEKPTHYR_                                         15.06326
    ## _AAAAALSGSPPQT[Phospho_STY_NL_98 (STY)]EKPTHYR_                                               NA
    ## _AAAAAVGPGLGS[Phospho_STY_NL_98 (STY)]GPGDS[Phospho_STY_NL_98 (STY)]PEGPEADAPER_        18.67455
    ##                                                                                  Exp2_vHC_CON_1
    ## _AAAAAAAAAAAAGDSDSWDADTFSMEDPVR_                                                       17.25051
    ## _AAAAALSGS[Phospho_STY_NL_98 (STY)]PPQTEKPTHYR_                                        15.46761
    ## _AAAAALSGSPPQT[Phospho_STY_NL_98 (STY)]EKPTHYR_                                        15.91756
    ## _AAAAAVGPGLGS[Phospho_STY_NL_98 (STY)]GPGDS[Phospho_STY_NL_98 (STY)]PEGPEADAPER_       18.66463
    ##                                                                                  Exp2_dHC_CON_1
    ## _AAAAAAAAAAAAGDSDSWDADTFSMEDPVR_                                                       17.36369
    ## _AAAAALSGS[Phospho_STY_NL_98 (STY)]PPQTEKPTHYR_                                        15.38014
    ## _AAAAALSGSPPQT[Phospho_STY_NL_98 (STY)]EKPTHYR_                                              NA
    ## _AAAAAVGPGLGS[Phospho_STY_NL_98 (STY)]GPGDS[Phospho_STY_NL_98 (STY)]PEGPEADAPER_       18.84291
    ##                                                                                  Exp2_vHC_SW30_1
    ## _AAAAAAAAAAAAGDSDSWDADTFSMEDPVR_                                                        17.33894
    ## _AAAAALSGS[Phospho_STY_NL_98 (STY)]PPQTEKPTHYR_                                         15.28279
    ## _AAAAALSGSPPQT[Phospho_STY_NL_98 (STY)]EKPTHYR_                                               NA
    ## _AAAAAVGPGLGS[Phospho_STY_NL_98 (STY)]GPGDS[Phospho_STY_NL_98 (STY)]PEGPEADAPER_        18.62069

lets access further information about each peptide (i.e phosphosites and probabilites)

    rowdat <- rowData(data_phos)
    head(rowdat)

    ## DataFrame with 6 rows and 7 columns
    ##                                                                                           Protein
    ##                                                                                          <factor>
    ## _AAAAAAAAAAAAGDSDSWDADTFSMEDPVR_                                                           Q66JS6
    ## _AAAAALSGS[Phospho_STY_NL_98 (STY)]PPQTEKPTHYR_                                            Q5FWH2
    ## _AAAAALSGSPPQT[Phospho_STY_NL_98 (STY)]EKPTHYR_                                            Q5FWH2
    ## _AAAAAVGPGLGS[Phospho_STY_NL_98 (STY)]GPGDS[Phospho_STY_NL_98 (STY)]PEGPEADAPER_           Q3UVL4
    ## _AAAC[Carbamidomethyl (C)]PS[Phospho_STY_NL_98 (STY)]S[Phospho_STY_NL_98 (STY)]PHKIPLSR_   Q812A2
    ## _AAAC[Carbamidomethyl (C)]PSS[Phospho_STY_NL_98 (STY)]PHKIPLSR_                            Q812A2
    ##                                                                                                StrippedSequence
    ##                                                                                                     <character>
    ## _AAAAAAAAAAAAGDSDSWDADTFSMEDPVR_                                                         AAAAAAAAAAAAGDSDSWDA..
    ## _AAAAALSGS[Phospho_STY_NL_98 (STY)]PPQTEKPTHYR_                                            AAAAALSGSPPQTEKPTHYR
    ## _AAAAALSGSPPQT[Phospho_STY_NL_98 (STY)]EKPTHYR_                                            AAAAALSGSPPQTEKPTHYR
    ## _AAAAAVGPGLGS[Phospho_STY_NL_98 (STY)]GPGDS[Phospho_STY_NL_98 (STY)]PEGPEADAPER_         AAAAAVGPGLGSGPGDSPEG..
    ## _AAAC[Carbamidomethyl (C)]PS[Phospho_STY_NL_98 (STY)]S[Phospho_STY_NL_98 (STY)]PHKIPLSR_        AAACPSSPHKIPLSR
    ## _AAAC[Carbamidomethyl (C)]PSS[Phospho_STY_NL_98 (STY)]PHKIPLSR_                                 AAACPSSPHKIPLSR
    ##                                                                                             Positions
    ##                                                                                           <character>
    ## _AAAAAAAAAAAAGDSDSWDADTFSMEDPVR_                                                                     
    ## _AAAAALSGS[Phospho_STY_NL_98 (STY)]PPQTEKPTHYR_                                          7;9;13;17;19
    ## _AAAAALSGSPPQT[Phospho_STY_NL_98 (STY)]EKPTHYR_                                          7;9;13;17;19
    ## _AAAAAVGPGLGS[Phospho_STY_NL_98 (STY)]GPGDS[Phospho_STY_NL_98 (STY)]PEGPEADAPER_                12;17
    ## _AAAC[Carbamidomethyl (C)]PS[Phospho_STY_NL_98 (STY)]S[Phospho_STY_NL_98 (STY)]PHKIPLSR_       6;7;14
    ## _AAAC[Carbamidomethyl (C)]PSS[Phospho_STY_NL_98 (STY)]PHKIPLSR_                                6;7;14
    ##                                                                                                  Probabilities
    ##                                                                                                    <character>
    ## _AAAAAAAAAAAAGDSDSWDADTFSMEDPVR_                                                                              
    ## _AAAAALSGS[Phospho_STY_NL_98 (STY)]PPQTEKPTHYR_                                                0.05;0.95;0;0;0
    ## _AAAAALSGSPPQT[Phospho_STY_NL_98 (STY)]EKPTHYR_                                          0.01;0.19;0.79;0.01;0
    ## _AAAAAVGPGLGS[Phospho_STY_NL_98 (STY)]GPGDS[Phospho_STY_NL_98 (STY)]PEGPEADAPER_                           1;1
    ## _AAAC[Carbamidomethyl (C)]PS[Phospho_STY_NL_98 (STY)]S[Phospho_STY_NL_98 (STY)]PHKIPLSR_                 1;1;0
    ## _AAAC[Carbamidomethyl (C)]PSS[Phospho_STY_NL_98 (STY)]PHKIPLSR_                                    0.01;0.99;0
    ##                                                                                          GeneSymbol
    ##                                                                                            <factor>
    ## _AAAAAAAAAAAAGDSDSWDADTFSMEDPVR_                                                             Eif3j2
    ## _AAAAALSGS[Phospho_STY_NL_98 (STY)]PPQTEKPTHYR_                                              Unkl  
    ## _AAAAALSGSPPQT[Phospho_STY_NL_98 (STY)]EKPTHYR_                                              Unkl  
    ## _AAAAAVGPGLGS[Phospho_STY_NL_98 (STY)]GPGDS[Phospho_STY_NL_98 (STY)]PEGPEADAPER_             Vps51 
    ## _AAAC[Carbamidomethyl (C)]PS[Phospho_STY_NL_98 (STY)]S[Phospho_STY_NL_98 (STY)]PHKIPLSR_     Srgap3
    ## _AAAC[Carbamidomethyl (C)]PSS[Phospho_STY_NL_98 (STY)]PHKIPLSR_                              Srgap3
    ##                                                                                          EarliestPos
    ##                                                                                            <integer>
    ## _AAAAAAAAAAAAGDSDSWDADTFSMEDPVR_                                                                   2
    ## _AAAAALSGS[Phospho_STY_NL_98 (STY)]PPQTEKPTHYR_                                                    7
    ## _AAAAALSGSPPQT[Phospho_STY_NL_98 (STY)]EKPTHYR_                                                    7
    ## _AAAAAVGPGLGS[Phospho_STY_NL_98 (STY)]GPGDS[Phospho_STY_NL_98 (STY)]PEGPEADAPER_                   2
    ## _AAAC[Carbamidomethyl (C)]PS[Phospho_STY_NL_98 (STY)]S[Phospho_STY_NL_98 (STY)]PHKIPLSR_         889
    ## _AAAC[Carbamidomethyl (C)]PSS[Phospho_STY_NL_98 (STY)]PHKIPLSR_                                  889
    ##                                                                                            PhosphoSites
    ##                                                                                             <character>
    ## _AAAAAAAAAAAAGDSDSWDADTFSMEDPVR_                                                                   None
    ## _AAAAALSGS[Phospho_STY_NL_98 (STY)]PPQTEKPTHYR_                                          13;15;19;23;25
    ## _AAAAALSGSPPQT[Phospho_STY_NL_98 (STY)]EKPTHYR_                                          13;15;19;23;25
    ## _AAAAAVGPGLGS[Phospho_STY_NL_98 (STY)]GPGDS[Phospho_STY_NL_98 (STY)]PEGPEADAPER_                  13;18
    ## _AAAC[Carbamidomethyl (C)]PS[Phospho_STY_NL_98 (STY)]S[Phospho_STY_NL_98 (STY)]PHKIPLSR_    894;895;902
    ## _AAAC[Carbamidomethyl (C)]PSS[Phospho_STY_NL_98 (STY)]PHKIPLSR_                             894;895;902

As you can see this now enables us to access information about
phosphosites and probabilities in this dataset
