Both 220308... and wtc11_full_gene... were taken from Evvie's directories, but should be recreated
In the screen design part of the pipeline


even_older_guide_targets.tsv comes from /oak/stanford/groups/engreitz/Users/jgalante/ENCODE_nonK562_Validation_Datasets/DC_TAP/CRISPRiScreen_w_Sceptre/resources/DC_TAP/old_guide_targets.tsv
where the naming convention comes from the fact that I used a similar "pre-merged" file from Evvie for the K562 DC TAP but "even older" is because in the K562 DC TAP I also used a post-merged old file


ChosenGenes.AllRegions.bed was copied from /oak/stanford/groups/engreitz/Users/dulguun/CRISPRDesign/220303_WTC11_Random_Screen_allTSS/01_ChooseRegions/ChosenGenes.AllRegions.bed

The ChosenGenes.DistalElements.bed was copied from /oak/stanford/groups/engreitz/Users/dulguun/CRISPRDesign/220310_WTC11_Random_Screen_controlDE_macs2/01_ChooseRegions/ChosenGenes.AllRegions.bed


# The guide design files were copied from [,] into the resources directory 
  # /oak/stanford/groups/engreitz/Users/dulguun/CRISPRDesign/220308_Random_Screen_FinalDesignFiles/220308_K562_Random_Screen_Crop.design.txt
  # /oak/stanford/groups/engreitz/Users/dulguun/CRISPRDesign/220308_Random_Screen_FinalDesignFiles/220308_WTC11_Random_Screen_Crop.design.txt
  
  
The sgOpti-DC file was modified by Dulguun to create the CROP file. Which is just a replacement in the "name" column so that sgOpti-DC became CROP. I'm currently using the older versino of the CROP file listd in the file path above which just replaces the name - not anything else i.e. the guide numbers are all the same