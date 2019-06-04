#' uCAREChemSuiteCLI
#' @title drug.class.deterministic
#' @description Takes structure data file (SDF) of candidate drug and predicts its drug class using deterministic model.
#' @param sdf input sdf file
#' @return predicted drug class of the candidate drug by deterministic model
#' @usage drug.class.deterministic("sdf")
#' @import ChemmineR
#' @import stats
#' @import utils
#' @import usethis
#' @examples{
#' example.class.deterministic<- system.file('extdata/example.sdf', package="uCAREChemSuiteCLI")
#' drug.class.deterministic(example.class.deterministic)
#' }
#' @export

drug.class.deterministic <- function(sdf)
{
  # Reading Query
  sdf_input <- read.SDFset(sdf[[1]])
  ac <- atomcountMA(read.SDFset(sdf[[1]]), addH=FALSE)

  # Acquisition of Attributes #############################################################
  # Extraction of Chemical content features
  carbon<-0
  hydrogen<-0
  oxygen<-0
  nitrogen<-0
  chlorine<-0
  sulfur<-0
  bromium<-0
  fluorine<-0
  ro<-0
  rings<-0
  aromatic_ring<-0
  quin<-0

  # Carbon Content
  if (("C" %in% colnames(ac)) == TRUE)
    carbon<- ac[,"C"]

  # Hydrogen Content
  if (("H" %in% colnames(ac)) == TRUE)
    hydrogen<- ac[,"H"]

  # Oxygen Content
  if (("O" %in% colnames(ac)) == TRUE)
    oxygen<- ac[,"O"]

  # Nitrogen Content
  if (("N" %in% colnames(ac)) == TRUE)
    nitrogen<- ac[,"N"]

  # Chlorine Content
  if (("Cl" %in% colnames(ac)) == TRUE)
    chlorine<- ac[,"Cl"]

  # Sulfur Content
  if (("S" %in% colnames(ac)) == TRUE)
    sulfur<- ac[,"S"]

  # Bromium Content
  if (("Br" %in% colnames(ac)) == TRUE)
    bromium<- ac[,"Br"]

  # Fluorine Content
  if (("F" %in% colnames(ac)) == TRUE)
    fluorine<- ac[,"F"]

  # Extraction of ring features
  ro <- rings(read.SDFset(sdf[[1]]),type="count",upper=1000, arom=TRUE)
  rings <-  rings(read.SDFset(sdf[[1]]), upper=Inf, type="all", arom=TRUE, inner=FALSE)

  # Quinolone condition starts here
  aromatic_rings<- list()
  aromatic_rings_final<- list()

  aliphatic_rings<- list()
  aliphatic_rings_final<- list()

  intersect_rings<- list()
  #intersect_rings_final<- list()

  if(!is.null(rings$AROMATIC))
  {
    for(i in 1:length(rings$AROMATIC))
    {
      if(rings$AROMATIC[i]==TRUE)
      {
        aromatic_rings<- rings$RINGS[i]
        aromatic_rings_final<- rbind(aromatic_rings_final, aromatic_rings)
      }
    }


    for(i in 1:length(rings$AROMATIC))
    {
      if(rings$AROMATIC[i]==FALSE)
      {
        aliphatic_rings<- rings$RINGS[i]
        aliphatic_rings_final<- rbind(aliphatic_rings_final, aliphatic_rings)
      }
    }

    if(as.numeric(length(rings$RINGS)) == as.numeric(length(aliphatic_rings_final)))
    {
      quin<-0
    }
    else
    {
      for(i in 1:length(aromatic_rings_final))
      {
        for(j in 1:length(aromatic_rings_final))
        {
          intersect_rings<-Reduce(intersect, list(aromatic_rings_final[[i]], aromatic_rings_final[[j]]))
          if(length(intersect_rings)==2)
          {
            quin<-1
          }
        }
      }
    }
  }else {quin <- 2}
  # Quinole condition ends here

  # Whether Aromatic ring is present or not
  if(!is.null(rings$AROMATIC))
  {
    for (i in 1: length(rings$AROMATIC))
    {
      if(rings$AROMATIC[[i]]=='TRUE' )
      {
        ring_p_or_a<-1
      }
    }
  }else{ring_p_or_a<-0}

  # $Acquisition of conditions for sulfonamide
  # Rule: Sulfur with two double bonded oxygen = TRUE or Not
  matrix<-  as.data.frame(conMA(read.SDFset(sdf[[1]]), exclude=c("H")))
  rowname_sul <- row.names(matrix)
  colname_sul <- colnames(matrix)
  row_sul_elements<-  data.frame(strsplit(rowname_sul,"_"))[1,]

  # Acquisition of Sulfur indices
  sul_index_S <- list()
  sul_index_current_S<- list()

  for(i in 1:length(row_sul_elements))
  {
    if(row_sul_elements[i]== 'S')
    {
      sul_index_current_S <- i
      sul_index_S <- rbind(sul_index_S,sul_index_current_S)
    }
  }

  # Acquisition of Oxygen indices
  sul_index_O <- list()
  sul_index_current_O<- list()

  for(i in 1:length(row_sul_elements))
  {
    if(row_sul_elements[i]== 'O')
    {
      sul_index_current_O <- i
      sul_index_O <- rbind(sul_index_O,sul_index_current_O)
    }
  }

  sul_matrix<-0
  if(sul_index_current_S[1] != "NULL")
  {
    output<- list()
    output_final<- list()
    for(i in 1: length(sul_index_S[,1]))
    {
      l<-0
      for(j in 1: length(sul_index_O[,1]))
      {
        a<- sul_index_S[,1][i]
        b<- sul_index_O[,1][j]
        if(matrix[a$sul_index_current_S[1],b$sul_index_current_O]==2)
        {
          l<-l+2
        }
      }
      output<- l
      output_final<- rbind(output_final,output)
    }

    for(i in 1: length(output_final))
    {
      if(output_final[i]==4)
      {
        sul_matrix<-4
      }
    }
  }

  # Rules for sulfonamide ends here
  ######################################## Formating needed after this #############################

  # Rule for Aminoglycoside (Whether glycosidic sugar present or not)
  aromatic_rings_final<- list()
  sugar.present<- 0

  if(!is.null(rings$AROMATIC))
  {
    for(i in 1:length(rings$AROMATIC))
    {
      if(rings$AROMATIC[i]==FALSE)
      {
        aromatic_rings<- rings$RINGS[i]
        ring.elements<-  data.frame(strsplit(aromatic_rings[[1]],"_"))[1,]


        aliphatic.oxygen.counter <- 0
        aliphatic.carbon.counter <- 0
        for(j in 1:length(ring.elements))
        {
          if("C" %in% ring.elements[,j] == TRUE)
          {
            aliphatic.carbon.counter <- aliphatic.carbon.counter + 1
          }
          else if("O" %in% ring.elements[,j] == TRUE)
          {
            aliphatic.oxygen.counter <- aliphatic.oxygen.counter + 1
          }
        }

        if(aliphatic.carbon.counter  == 5 && aliphatic.oxygen.counter == 1 && length(ring.elements)==6)
        {
          sugar.present <- 1
        }
      }
    }
  }


  #Functional groups common in Aminoglycosides viz. ROR and ROH
  functional.group<- groups(read.SDFset(sdf[[1]]), groups="fctgroup", type="countMA")

  imp.func.group.aminoglycoside <- 0
  if(functional.group["ROH"] != 0 && functional.group["ROR"] !=0)
  {
    imp.func.group.aminoglycoside <- 1
  }




  #Rule for Nitrofurantoin
  #Acquisition of Carbon number in the aromatic ring of nitrofurans

  if(!is.null(rings$AROMATIC))
  {
    top1<- list()
    top1_0<- list()
    for (i in 1:length(rings$AROMATIC))
    {
      if(rings$AROMATIC[i]== TRUE)
      {
        X <- list(rings$RINGS[i])[[1]][[1]]

        if(length(X)==5)
        {
          Y <- data.frame(strsplit(X,"_"))[1,]
          Z <- rowSums(Y == "C")
          J <- rowSums(Y== "O")

          output_O<- J[[1]]
          output_C<- Z[[1]]

          top<- as.data.frame(output_C)
          top_O<- as.data.frame(output_O)

          top1<- rbind(top1,top)
          top1_0<- rbind(top1_0, top_O)
        }
      }
      k<-0
      if(!is.null(ncol(top1)))
      {
        for(j in 1: ncol(top1))
        {
          if(top1[j]== 4 && top1_0[j]==1)
          {
            k<-4
          }
        }
      }
      else
      {
        k<-0
      }

      ring_nitrofuran_C<-k
    }
  }else(ring_nitrofuran_C <- 0)

  #Acquisition of Carbon content in aliphatic rings
  if(!is.null(rings$AROMATIC) && !is.null(rings$RINGS))
  {
    top1<- list()
    for (i in 1:length(rings$AROMATIC))
    {
      if(rings$AROMATIC[i]== FALSE)
      {
        X <- list(rings$RINGS[i])[[1]][[1]]
        Y <- data.frame(strsplit(X,"_"))[1,]
        Z <- rowSums(Y == "C")

        output_C<- Z[[1]]
        top<- as.data.frame(output_C)
        top1<- rbind(top1,top)
      }
    }
    if(!is.null(nrow(top1)))
    {
      ring_C<-0
      for(j in 1: nrow(top1))
      {
        if(top1[[1]][j]== 3)
        {
          ring_C<-3
        }
      }
    }else
    { ring_C <- 0}

  } else
  { ring_C <- 0}

  #Acquisition of Nitrogen content in aliphatic rings
  if(!is.null(rings$AROMATIC) && !is.null(rings$RINGS))
  {
    top_N<- list()
    top1_N<- list()
    for (i in 1:length(rings$AROMATIC))
    {
      if(rings$AROMATIC[i] == FALSE)
      {
        X <- list(rings$RINGS[i])[[1]][[1]]
        Y <- data.frame(strsplit(X,"_"))[1,]
        Z <- rowSums(Y == "N")

        output_N<- Z[[1]]
        top_N<- as.data.frame(output_N)
        top1_N<- rbind(top1_N,top_N)
      }
    }
    ring_N<-0

    if(!is.null(nrow(top1_N)))
    {
      for(j in 1: nrow(top1_N))
      {
        if(top1_N[[1]][j]== 1)
        {
          ring_N<-1
        }
      }
    }else {ring_N<-0}


  }else{ ring_N<-0 }

  #Total Rings
  if(!is.null(ro[1]))
  {
    total_ring<- ro[1]
  }

  #Aromatic rings
  if(!is.null(ro[2]))
  {
    aromatic_ring<-ro[2]
  }

  # Acquisition of Attributes ends here #############################################################


  # Prediction (Deterministic) starts here ######
  # Rule for Pyridopyrimidine
  if(carbon == 14 && hydrogen == 17 &&  nitrogen == 5  && oxygen ==3){
    predicted_class<- "Pyridopyrimidine"
  }

  # Rule for Aminoquinone
  else if(carbon == 25 && hydrogen == 22 &&  nitrogen == 4  && oxygen ==8){
    predicted_class<- "Aminoquinone"
  }

  #Rule for Acriflavine
  else if(carbon == 27 && chlorine ==1 && hydrogen == 25 &&  nitrogen == 6 && total_ring==12 && aromatic_ring ==12){
    predicted_class<- "Acriflavine"
  }

  #Rule for Aminocoumarin
  else if(carbon ==31 && hydrogen == 36 && nitrogen==2 && oxygen ==11 && total_ring== 5 && aromatic_ring==4)
  {
    predicted_class<- "Aminocoumarin"
  }

  #Rule for Anisole
  else if (carbon == 14 && hydrogen == 18 && nitrogen == 4 && oxygen == 3 && total_ring == aromatic_ring)
  {
    predicted_class<- "Anisole"
  }

  #Rule for Anthracycline
  else if(carbon == 27 && hydrogen >= 29 && oxygen ==10  && nitrogen ==1 && chlorine ==0 && total_ring ==11&& aromatic_ring >=2 ){
    predicted_class<- "Anthracycline"
  }

  #Rule for Benzalkonium
  else if(bromium ==1&& carbon == 21 && hydrogen  >= 38 && oxygen ==0  && nitrogen ==1 && chlorine ==0 && total_ring ==1 && aromatic_ring ==1 ){
    predicted_class<-"Benzalkonium"
  }

  #Rule for Chloramphenicol
  else if(carbon == 11 && hydrogen == 12 && oxygen ==5  && nitrogen ==2 && chlorine ==2 && total_ring ==1 && aromatic_ring ==1 ){
    predicted_class<- "Chloramphenicol"
  }

  #Rule for Florfenicol
  else if(carbon != 0 && hydrogen != 0 && chlorine !=0 && fluorine !=0 && nitrogen !=0 && oxygen !=0 && sulfur !=0    && total_ring ==1 && aromatic_ring ==1 ){
    predicted_class<-"Florfenicol"
  }

  #Rule for Rhodamine
  else if(carbon == 28 && hydrogen == 31 && chlorine ==1 && fluorine ==0 && nitrogen ==2 && oxygen ==3 && sulfur ==0    && total_ring ==7 && aromatic_ring ==7 ){
    predicted_class<- "Rhodamine"
  }

  #Rule for Thiolactomycin
  else if(carbon == 11 && hydrogen == 14 && chlorine ==0 && fluorine ==0 && nitrogen ==0 && oxygen ==2 && sulfur ==1    && total_ring ==1 && aromatic_ring ==0 ){
    predicted_class<- "Thiolactomycin"
  }

  #Rule for Beta Lactam
  else if(carbon >= 12 && carbon <=23 && hydrogen >= 16 && hydrogen <= 27 && oxygen >=4 && oxygen <=9 && nitrogen >=2 && nitrogen <=8 && sulfur >=1 && sulfur <=3 && ring_C ==3 && ring_N ==1){
    predicted_class<- "Beta Lactam"
  }

  #Rule for Sulfonamide
  else if(carbon != 0 && hydrogen != 0 && oxygen != 0 && nitrogen !=0 && sulfur !=0 && sul_matrix == 4 && ring_p_or_a == 1){
    predicted_class<- "Sulfonamide"
  }


  #Rule for Fluoroquinolone
  else if(carbon != 0 && hydrogen != 0 && oxygen == 3 && nitrogen ==3 && fluorine !=0 && quin==1){
    predicted_class<- "Fluoroquinolone"
  }

  #Rule for Polyketide_Erythromycin
  else if(carbon == 37 && hydrogen == 67 && nitrogen ==1 && oxygen == 13 && fluorine ==0){
    predicted_class<- "Polyketide"
  }

  #Rule for Polyketide_Clarithromycin
  else if(carbon == 38 && hydrogen == 69 && nitrogen ==1 && oxygen == 13 && fluorine ==0){
    predicted_class<- "Polyketide"
  }

  #Rule for Polyketide_Cethromycin
  else if(carbon == 42 && hydrogen == 59 && nitrogen ==3 && oxygen == 10 && fluorine ==0){
    predicted_class<- "Polyketide"
  }

  #Rule for Polyketide_Telithromycin
  else if(carbon == 43 && hydrogen == 65 && nitrogen ==5 && oxygen == 10 && fluorine ==0){
    predicted_class<- "Polyketide"
  }

  #Rule for Polyketide_Rifampin
  else if(carbon == 43 && hydrogen == 58 && nitrogen ==4 && oxygen == 12 && fluorine ==0){
    predicted_class<- "Polyketide"
  }

  #Rule for Polyketide_Tetracycline
  else if(carbon == 22 && hydrogen == 24 && nitrogen == 2 && oxygen == 8 && fluorine ==0){
    predicted_class<- "Polyketide"
  }

  #Rule for Polyketide_Minocycline
  else if(carbon == 23 && hydrogen == 27 && nitrogen == 3 && oxygen == 7 && fluorine ==0){
    predicted_class<- "Polyketide"
  }

  #Rule for Quinolone
  else if(carbon != 0 && hydrogen != 0 && oxygen == 3 && nitrogen ==2 && fluorine ==0 && quin==1){
    predicted_class<- "Quinolone"
  }

  #Rule for Peptide drug Bicyclomycin (C12H18N2O7)
  else if(carbon ==12 && hydrogen ==18 && nitrogen ==2 && oxygen ==7){
    predicted_class<- "Peptide drug"
  }

  #Rule for Peptide drug Polymyxin B (C56H100N16O17S)
  else if(carbon ==56 && hydrogen ==98 && nitrogen ==16 && oxygen ==13){
    predicted_class<- "Peptide drug"
  }

  #Rule for Aminoglycoside
  else if(carbon != 0 && hydrogen >= 2 * oxygen && nitrogen !=0 && chlorine ==0 && sugar.present ==1 && imp.func.group.aminoglycoside == 1){
    predicted_class<- "Aminoglycosides"
  }

  #Rule for Nitrofurans
  else if(nitrogen !=0 && ring_nitrofuran_C == 4){
    predicted_class<- "Nitrofurans"
  }

  else
  {
    predicted_class<- "The drug is not found in the list! Please use Stochastic model"

  }

  return(predicted_class)


}
# Deterministic model ends here



# Stochastic class predition starts here
#' uCAREChemSuiteCLI
#' @title drug.class.stochastic
#' @description Takes structure data file (SDF) of candidate drug, Nearest Neighbor value and threshold similarity score to predict its drug class using stochastic model.
#' @param sdf input sdf file
#' @param NearestNeighbor Nearest Neighbor = 1, 3
#' @param Threshold Threshold = 0.25, 0.3, 0.35, 0.4
#' @return Predicted drug class of the candidate drug using Nearest Neighbor algorithm
#' @usage drug.class.stochastic("sdf", "NearestNeighbor", "Threshold")
#' @import ChemmineR
#' @import stats
#' @import utils
#' @import usethis
#' @examples{
#' example.class.stochastic<- system.file('extdata/example.sdf', package="uCAREChemSuiteCLI")
#' drug.class.stochastic(example.class.stochastic,"3","0.25")
#' }
#' @export

drug.class.stochastic <- function(sdf, NearestNeighbor, Threshold){


  # Reading the Database
  db_antibiotics<- system.file('extdata/all_sdf_names.sdf', package="uCAREChemSuiteCLI")
  antibiotics <- read.SDFset(db_antibiotics)
  cid(antibiotics) <-  makeUnique(sdfid(antibiotics))
  apset <- sdf2ap(antibiotics)
  Score<- 0
  Drug_Name<-0


  # Reading Query
  sdf_input <- read.SDFset(sdf[[1]])
  sdf_input_ap<- sdf2ap(sdf_input)

  df<-cmp.search(apset, sdf_input_ap, type=3, cutoff = 77, quiet=TRUE)
  colnames(df)<- c("Index","Drug_Name","Score")

  # Clustering drug against database
  subset_of_drugs_with_equal_or_less_than_score <- subset(df, Score >= Threshold)
  drugs_falling_under_threshold<-head(subset_of_drugs_with_equal_or_less_than_score, n=as.numeric(NearestNeighbor))

  if(nrow(subset_of_drugs_with_equal_or_less_than_score) != 0)
  {
    # Finding the class of drug using Stochastic model
    antibiotic_dataset<- system.file('extdata/Antibiotic_data_set.csv', package="uCAREChemSuiteCLI")
    db<- read.csv2(antibiotic_dataset ,header = TRUE, sep=",")
    drug_list<- drugs_falling_under_threshold[2]

    if (as.numeric(table(subset_of_drugs_with_equal_or_less_than_score["Score"] >=Threshold)["TRUE"]) < as.numeric(NearestNeighbor))
    {
      subset_of_drug_class<- data.frame()
      subset_of_drug_class_final<- data.frame()

      for(i in 1:length(drug_list[,1]))
      {
        for(j in 1:length(db[,1]))
        {
          drug<- drug_list[[1]][i]
          db_drug<- db[[1]][j]


          db_drug_final <- gsub(" ", "", db_drug[1])
          drug_final <- gsub("_", "", drug[1])

          if(db_drug_final[1] %in% drug_final[1])
          {
            subset_of_drug_class<-db[j,]
            subset_of_drug_class_final<- rbind(subset_of_drug_class_final,subset_of_drug_class)
          }
        }
      }

      #Concatanated dataframe of drugname, drug class and similarity socre
      concat_table_drug_Class_score<- cbind(unique(subset_of_drug_class_final[,1:2]),drugs_falling_under_threshold[,3])

      #Prediction of drug class
      if(as.numeric(NearestNeighbor) == 1)
      {
        predicted_class<- c("Drug class:",as.character(concat_table_drug_Class_score[1,2]))
      }
      else if(as.numeric(NearestNeighbor) == 3)
      {
        if(as.numeric(concat_table_drug_Class_score[1,3])==1)
        {
          predicted_class<- c("Drug class:",as.character(concat_table_drug_Class_score[1,2]))
        }
        else
        {
          drug_class_frequency<- aggregate(data.frame(count = concat_table_drug_Class_score[,2]), list(value = concat_table_drug_Class_score[,2]), length)
          if(as.numeric(drug_class_frequency[1,2])>=2)
          {
            predicted_class<- c("Drug class:",as.character(drug_class_frequency[1,1]))
          }
          else
          {
            predicted_class<- c("Drug class:",as.character(drug_class_frequency[1,1]))
          }
        }
      }
    }
    else
    {
      neighbor <-subset_of_drugs_with_equal_or_less_than_score[1,2]
      fetchted_row<-   subset(db, Drug_Name %in% as.character(neighbor))
      predicted_class<-c("Drug class:", as.character(fetchted_row["Drug_class"][1,1]))
    }
  }
  else
  {
    predicted_class<- "Not enough neighbors"
  }
  #predicted_class<- db_drug[1]
  return(predicted_class)
  # Stochastic model Ends here ###############################################
}

# Resistome prediction starts here
#' uCAREChemSuiteCLI
#' @title drug.resistome.deterministic
#' @description Takes structure data file (SDF) of candidate drug to predicts its resistome using deterministic model.
#' @param sdf input sdf file
#' @param Organism Escherichia coli = 1, Pseudomonas aeruginosa = 2
#' @return Predicted resistome of the candidate drug using deterministic model
#' @usage drug.resistome.deterministic("sdf", "Organism")
#' @import ChemmineR
#' @import stats
#' @import utils
#' @import usethis
#' @examples{
#' example.resistome.deterministic<- system.file('extdata/example.sdf', package="uCAREChemSuiteCLI")
#' drug.resistome.deterministic(example.resistome.deterministic, "1")
#' }
#' @export


drug.resistome.deterministic <- function(sdf, Organism)
{

  drug_predicted_class<- drug.class.deterministic(sdf)

  # Reading the Database
  if(Organism == 1)
  {
    antibiotic_dataset<- system.file('extdata/Antibiotic_data_set_Ecoli.csv', package="uCAREChemSuiteCLI")
    db<- read.csv2(antibiotic_dataset ,header = TRUE, sep=",")
    db_drug_class<-db["Drug_class"]
  } else if(Organism == 2)
  {antibiotic_dataset<- system.file('extdata/Antibiotic_data_set_Paeruginosa.csv', package="uCAREChemSuiteCLI")
  db<- read.csv2(antibiotic_dataset ,header = TRUE, sep=",")
  db_drug_class<-db["Drug_class"]
  }

  drug_resistome<- list()

  for(i in 1: nrow(db_drug_class))
  {
    if(db_drug_class[[1]][i] %in% drug_predicted_class)
    {
      drug_resistome_raw <- db[i,]
      drug_resistome<- rbind(drug_resistome, drug_resistome_raw)
    }

  }

  return(drug_resistome)


}


#' uCAREChemSuiteCLI
#' @title drug.resistome.stochastic
#' @description Takes structure data file (SDF) of candidate drug to predict its resistome using stochastic model.
#' @param sdf input sdf file
#' @param NearestNeighbor Nearest Neighbor = 1, 3
#' @param Threshold Threshold = 0.25, 0.3, 0.35, 0.4
#' @param Organism Escherichia coli = 1, Pseudomonas aeruginosa = 2
#' @return Predicted resistome of the candidate drug using Nearest Neighbor algorithm
#' @usage drug.resistome.stochastic("sdf", "NearestNeighbor", "Threshold", "Organism")
#' @import ChemmineR
#' @import stats
#' @import utils
#' @import usethis
#' @examples{
#' example.resistome.stochastic<- system.file('extdata/example.sdf', package="uCAREChemSuiteCLI")
#' drug.resistome.stochastic(example.resistome.stochastic, "3", "0.25", "1")
#' }
#' @export

drug.resistome.stochastic <- function(sdf, NearestNeighbor, Threshold, Organism)
{

  drug_predicted_class<- as.data.frame(drug.class.stochastic(sdf, NearestNeighbor, Threshold))
  if(is.na(drug_predicted_class[[1]][2]))
  {
    drug_resistome<- "Not enough neighbors"
  }
  else
  {
    # Reading the Database
    if(Organism == 1){
      antibiotic_dataset<- system.file('extdata/Antibiotic_data_set_Ecoli.csv', package="uCAREChemSuiteCLI")
      db<- read.csv2(antibiotic_dataset ,header = TRUE, sep=",")
      db_drug_class<-db["Drug_class"]
    }else if(Organism == 2){
      antibiotic_dataset<- system.file('extdata/Antibiotic_data_set_Paeruginosa.csv', package="uCAREChemSuiteCLI")
      db<- read.csv2(antibiotic_dataset ,header = TRUE, sep=",")
      db_drug_class<-db["Drug_class"]
    }
    drug_resistome<- list()

    for(i in 1: nrow(db_drug_class))
    {
      if(db_drug_class[[1]][i] %in% drug_predicted_class[[1]][2])
      {
        drug_resistome_raw <- db[i,]
        drug_resistome<- rbind(drug_resistome, drug_resistome_raw)
      }
    }
  }

  return(drug_resistome)


}



