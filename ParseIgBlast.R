#This function goes through the IgBlast result file line by line and changes the results into tabular format for downstream analysis

#Input - IgBlast result file

#Output - A file containing tabular for of the results

ParseIgBlast <- function(file_name) {
  # create a folder to store the results
  if(!file.exists("Parse_IgBlast"))
  { dir.create("Parse_IgBlast")}
  #read the IgBlast results
  rawHTML <- paste(readLines(file_name))
  #calculate number of lines
  long <-length(rawHTML)

  #count the number of sequences
  number_of_sequences = 0
  for (b in 1:long)
  {
    condition <- grep("Effective",rawHTML[b] ) #Effective is present once for each query sequence result
    if(length(condition) == 1 )
        {
            number_of_sequences = number_of_sequences + 1
        }
  }
  rm(b, condition)
  
  #create a table to store information after parsing - one for each - heavy and light chain
  All_data_Heavy <- matrix( , nrow = number_of_sequences, ncol = 68)
  colnames(All_data_Heavy) = c("Query","Length","VH_top_match","Score_VH","E_value_VH","DH_top_match",
                               "Score_DH","E_value_DH","JH_top_match","Score_JH","E_value_JH",
                               "Chain_type","Stop_codon","V-Jframe","Productive","Strand","V_region_end",
                               "V-D_junction","D_region","D-J_junction","J-region_start","V_FR1-IMGT_from","V_FR1-IMGT_to","V_FR1-IMGT_length",
                               "V_FR1-IMGT_match","V_FR1-IMGT_mismatch","V_FR1-IMGT_gaps","V_FR1-IMGT_identity_percent",
                               "CDR1-IMGT_from","CDR1-IMGT_to","CDR1-IMGT_length","CDR1-IMGT_match","CDR1-IMGT_mismatch",
                               "CDR1-IMGT_gaps","CDR1-IMGT_identity_percent","FR2-IMGT_from","FR2-IMGT_to","FR2-IMGT_length",
                               "FR2-IMGT_match","FR2-IMGT_mismatch","FR2-IMGT_gaps","FR2-IMGT_identity_percent","CDR2-IMGT_from",
                               "CDR2-IMGT_to","CDR2-IMGT_length","CDR2-IMGT_match","CDR2-IMGT_mismatch","CDR2-IMGT_gaps",
                               "CDR2-IMGT_identity_percent","FR3-IMGT_from","FR3-IMGT_to","FR3-IMGT_length","FR3-IMGT_match",
                               "FR3-IMGT_mismatch","FR3-IMGT_gaps","FR3-IMGT_identity_percent","CDR3-IMGT_germline_from",
                               "CDR3-IMGT_germline_to","CDR3-IMGT_germline_length","CDR3-IMGT_germline_match",
                               "CDR3-IMGT_germline_mismatch","CDR3-IMGT_germline_gaps","CDR3-IMGT_germline_identity_percent",
                               "Total-IMGT_length","Total-IMGT_match","Total-IMGT_mismatch","Total-IMGT_gaps",
                               "Total-IMGT_identity_percent")
  All_data_Light <- matrix( , nrow = number_of_sequences, ncol = 65)
  colnames(All_data_Light) = c("Query","Length","VL_top_match","Score_VL","E_value_VL","JL_top_match","Score_JL","E_value_JL",
                               "Chain_type","Stop_codon","V-Jframe","Productive","Strand","V_region_end",
                               "V-J_junction","J-region_start","CDR3_nucleotide_seq",
                               "CDR3_aminoacid_seq","V_FR1-IMGT_from","V_FR1-IMGT_to","V_FR1-IMGT_length",
                               "V_FR1-IMGT_match","V_FR1-IMGT_mismatch","V_FR1-IMGT_gaps","V_FR1-IMGT_identity_percent",
                               "CDR1-IMGT_from","CDR1-IMGT_to","CDR1-IMGT_length","CDR1-IMGT_match","CDR1-IMGT_mismatch",
                               "CDR1-IMGT_gaps","CDR1-IMGT_identity_percent","FR2-IMGT_from","FR2-IMGT_to","FR2-IMGT_length",
                               "FR2-IMGT_match","FR2-IMGT_mismatch","FR2-IMGT_gaps","FR2-IMGT_identity_percent","CDR2-IMGT_from",
                               "CDR2-IMGT_to","CDR2-IMGT_length","CDR2-IMGT_match","CDR2-IMGT_mismatch","CDR2-IMGT_gaps",
                               "CDR2-IMGT_identity_percent","FR3-IMGT_from","FR3-IMGT_to","FR3-IMGT_length","FR3-IMGT_match",
                               "FR3-IMGT_mismatch","FR3-IMGT_gaps","FR3-IMGT_identity_percent","CDR3-IMGT_germline_from",
                               "CDR3-IMGT_germline_to","CDR3-IMGT_germline_length","CDR3-IMGT_germline_match",
                               "CDR3-IMGT_germline_mismatch","CDR3-IMGT_germline_gaps","CDR3-IMGT_germline_identity_percent",
                               "Total-IMGT_length","Total-IMGT_match","Total-IMGT_mismatch","Total-IMGT_gaps",
                               "Total-IMGT_identity_percent")
  
  main_file_line_count = 1 #Keeps a tab of number of lines parsed in the original results file
  seq_number = 0   #Keeps tab on number of sequences parsed
  heavy_num = 0   #Keeps tab on number of heavy chain sequences parsed
  light_num = 0    #Keeps tab on number of light chain sequences parsed
  
  #This loope runs till the last line of the document or the total number of sequences
  while (main_file_line_count != long && seq_number != number_of_sequences)
  {
    
        #look for start of result of each query
        condition <- grep("Query=",rawHTML[main_file_line_count] )  #the term "Query=" is present at the begining of the IgBlast results
        seq_number = seq_number + 1   #count for sequence number processed
        while (length(condition) !=1)
        {
                main_file_line_count = main_file_line_count + 1
                condition <- grep("Query=",rawHTML[main_file_line_count])
        }
        
        #look for end of result of each query
        condition <- grep("Effective",rawHTML[main_file_line_count] )  #the term "Effective" is present at the end of the IgBlast results for each sequence
        singleton = character(length = 0)  # stores results of one sequence
        while (length(condition) !=1)
        {
            condition <- grep("Effective",rawHTML[main_file_line_count] )
            singleton = rbind(singleton,rawHTML[main_file_line_count])
            main_file_line_count = main_file_line_count + 1
        }
        
        line_count_singleton <- length(singleton) # count number of lines in results of singleton sequence
        c = 1   #count to read line one by one
    
        # check whether its a "Heavy" or "Kappa" or "Lambda" or "No hit"
        while(!((grepl("lcl|IG", singleton[c]))||(grepl("No hits found", singleton[c])))) {c = c + 1 }
        if(grepl("IGHV", singleton[c]))
        {
            Chain <- ("Heavy")
        } else if(grepl("IGKV", singleton[c]))
        {
            Chain <- ("Kappa")
        } else if(grepl("IGLV", singleton[c]))
        {
            Chain <- ("Lambda")
        } else if(grepl("No hits found", singleton[c]))
        {
            Chain <- ("no hit")
            rm(c, line_count_singleton, singleton)
        }
    
        #look for markers in each line and store corresponding results into the tables
        if(Chain == "Heavy"||Chain == "Kappa"|| Chain == "Lambda")
        {
            v_string = strsplit(singleton[c], "\\s+")
            
            #parse and save results only if igBlast score is greater than 40
            if(as.numeric(v_string[[1]][2]) > 40)
            {
                #if it is a heavy chain
                if(Chain == "Heavy")
                {
                    heavy_num = heavy_num + 1
                    All_data_Heavy[heavy_num,"Query"] = gsub("Query= ","",singleton[1])
                    All_data_Heavy[heavy_num,"Length"] = gsub("Length=","", singleton[3])
                    v_string = strsplit(singleton[c], "\\s+")
                    All_data_Heavy[heavy_num,"Score_VH"] = v_string[[1]][2]
                    All_data_Heavy[heavy_num,"E_value_VH"] = v_string[[1]][3]
                    while(grepl("IGHV", singleton[c])||grepl("IGKV", singleton[c])||grepl("IGLV", singleton[c])) {c = c + 1}
                    if(grepl("IGHD", singleton[c]))
                    {
                        v_string = strsplit(singleton[c], "\\s+")
                        All_data_Heavy[heavy_num,"Score_DH"] = v_string[[1]][2]
                        All_data_Heavy[heavy_num,"E_value_DH"] = v_string[[1]][3]
                        while(grepl("IGHD", singleton[c])) {c = c + 1}
                    }
                    if(grepl("IGHJ", singleton[c])||grepl("IGKJ", singleton[c])||grepl("IGLJ", singleton[c]))
                    {
                        v_string = strsplit(singleton[c], "\\s+")
                        All_data_Heavy[heavy_num,"Score_JH"] = v_string[[1]][2]
                        All_data_Heavy[heavy_num,"E_value_JH"] = v_string[[1]][3]
                        while(grepl("IGHJ", singleton[c])||grepl("IGKJ", singleton[c])||grepl("IGLJ", singleton[c])) {c = c + 1}
                    }
          
                    while(!(grepl("Chain\ type", singleton[c]))){c = c + 1 }
                    v_string = strsplit(singleton[c+1], "\t")
                    All_data_Heavy[heavy_num,"VH_top_match"] = v_string[[1]][1]
                    All_data_Heavy[heavy_num,"DH_top_match"] = v_string[[1]][2]
                    All_data_Heavy[heavy_num,"JH_top_match"] = v_string[[1]][3]
                    All_data_Heavy[heavy_num,"Chain_type"]= v_string[[1]][4]
                    All_data_Heavy[heavy_num,"Stop_codon"]= v_string[[1]][5]
                    All_data_Heavy[heavy_num,"V-Jframe"]= v_string[[1]][6]
                    All_data_Heavy[heavy_num,"Productive"]= v_string[[1]][7]
                    All_data_Heavy[heavy_num,"Strand"]= v_string[[1]][8]
          
                    while(!(grepl("junction\ details", singleton[c]))){c = c + 1 }
                    v_string = strsplit(singleton[c + 1], "\t")
                    All_data_Heavy[heavy_num,"V_region_end"]= v_string[[1]][1]
                    All_data_Heavy[heavy_num,"V-D_junction"]= v_string[[1]][2]
                    All_data_Heavy[heavy_num,"D_region"]=v_string[[1]][3]
                    All_data_Heavy[heavy_num,"D-J_junction"]= v_string[[1]][4]
                    All_data_Heavy[heavy_num,"J-region_start"]= v_string[[1]][5]
                    singleton[c]
                    while(c <= line_count_singleton)
                    {
                        c = c + 3
                        if(grepl("Alignment\ summary",singleton[c]))
                        {
                            c = c + 1
                            if(grepl("FR1-IMGT",singleton[c]))
                            {
                                v_string = strsplit(singleton[c], "\t")
                                All_data_Heavy[heavy_num,"V_FR1-IMGT_from"]= v_string[[1]][2]
                                All_data_Heavy[heavy_num,"V_FR1-IMGT_to"]= v_string[[1]][3]
                                All_data_Heavy[heavy_num,"V_FR1-IMGT_length"]= v_string[[1]][4]
                                All_data_Heavy[heavy_num,"V_FR1-IMGT_match"]= v_string[[1]][5]
                                All_data_Heavy[heavy_num,"V_FR1-IMGT_mismatch"]=v_string[[1]][6]
                                All_data_Heavy[heavy_num,"V_FR1-IMGT_gaps"]= v_string[[1]][7]
                                All_data_Heavy[heavy_num,"V_FR1-IMGT_identity_percent"]=v_string[[1]][8]
                                c= c + 1
                            }
                            if(grepl("CDR1-IMGT",singleton[c]))
                            {
                                v_string = strsplit(singleton[c], "\t")
                                All_data_Heavy[heavy_num,"CDR1-IMGT_from"]= v_string[[1]][2]
                                All_data_Heavy[heavy_num,"CDR1-IMGT_to"]= v_string[[1]][3]
                                All_data_Heavy[heavy_num,"CDR1-IMGT_length"]= v_string[[1]][4]
                                All_data_Heavy[heavy_num,"CDR1-IMGT_match"]= v_string[[1]][5]
                                All_data_Heavy[heavy_num,"CDR1-IMGT_mismatch"]= v_string[[1]][6]
                                All_data_Heavy[heavy_num,"CDR1-IMGT_gaps"]= v_string[[1]][7]
                                All_data_Heavy[heavy_num,"CDR1-IMGT_identity_percent"]= v_string[[1]][8]
                                c= c + 1
                            }
                            if(grepl("FR2-IMGT",singleton[c]))
                            {
                                v_string = strsplit(singleton[c], "\t")
                                All_data_Heavy[heavy_num,"FR2-IMGT_from"]= v_string[[1]][2]
                                All_data_Heavy[heavy_num,"FR2-IMGT_to"]= v_string[[1]][3]
                                All_data_Heavy[heavy_num,"FR2-IMGT_length"]= v_string[[1]][4]
                                All_data_Heavy[heavy_num,"FR2-IMGT_match"]= v_string[[1]][5]
                                All_data_Heavy[heavy_num,"FR2-IMGT_mismatch"]= v_string[[1]][6]
                                All_data_Heavy[heavy_num,"FR2-IMGT_gaps"]= v_string[[1]][7]
                                All_data_Heavy[heavy_num,"FR2-IMGT_identity_percent"]= v_string[[1]][8]
                                c= c + 1
                            }
                            if(grepl("CDR2-IMGT",singleton[c]))
                            {
                                v_string = strsplit(singleton[c], "\t")
                                All_data_Heavy[heavy_num,"CDR2-IMGT_from"]= v_string[[1]][2]
                                All_data_Heavy[heavy_num,"CDR2-IMGT_to"]= v_string[[1]][3]
                                All_data_Heavy[heavy_num,"CDR2-IMGT_length"]= v_string[[1]][4]
                                All_data_Heavy[heavy_num,"CDR2-IMGT_match"]= v_string[[1]][5]
                                All_data_Heavy[heavy_num,"CDR2-IMGT_mismatch"]= v_string[[1]][6]
                                All_data_Heavy[heavy_num,"CDR2-IMGT_gaps"]= v_string[[1]][7]
                                All_data_Heavy[heavy_num,"CDR2-IMGT_identity_percent"]= v_string[[1]][8]
                                c= c + 1
                            }
                            if(grepl("FR3-IMGT",singleton[c]))
                            {
                                v_string = strsplit(singleton[c], "\t")
                                All_data_Heavy[heavy_num,"FR3-IMGT_from"]= v_string[[1]][2]
                                All_data_Heavy[heavy_num,"FR3-IMGT_to"]= v_string[[1]][3]
                                All_data_Heavy[heavy_num,"FR3-IMGT_length"]= v_string[[1]][4]
                                All_data_Heavy[heavy_num,"FR3-IMGT_match"]= v_string[[1]][5]
                                All_data_Heavy[heavy_num,"FR3-IMGT_mismatch"]= v_string[[1]][6]
                                All_data_Heavy[heavy_num,"FR3-IMGT_gaps"]= v_string[[1]][7]
                                All_data_Heavy[heavy_num,"FR3-IMGT_identity_percent"]= v_string[[1]][8]
                                c= c + 1
                            }
                            if(grepl("CDR3-IMGT",singleton[c]))
                            {
                                v_string = strsplit(singleton[c], "\t")
                                All_data_Heavy[heavy_num,"CDR3-IMGT_germline_from"]= v_string[[1]][2]
                                All_data_Heavy[heavy_num,"CDR3-IMGT_germline_to"]= v_string[[1]][3]
                                All_data_Heavy[heavy_num,"CDR3-IMGT_germline_length"]= v_string[[1]][4]
                                All_data_Heavy[heavy_num,"CDR3-IMGT_germline_match"]= v_string[[1]][5]
                                All_data_Heavy[heavy_num,"CDR3-IMGT_germline_mismatch"]= v_string[[1]][6]
                                All_data_Heavy[heavy_num,"CDR3-IMGT_germline_gaps"]= v_string[[1]][7]
                                All_data_Heavy[heavy_num,"CDR3-IMGT_germline_identity_percent"]= v_string[[1]][8]
                                c= c + 1
                            }
                            if(grepl("Total",singleton[c]))
                            {
                                v_string = strsplit(singleton[c], "\t")
                                All_data_Heavy[heavy_num,"Total-IMGT_length"]= v_string[[1]][4]
                                All_data_Heavy[heavy_num,"Total-IMGT_match"]= v_string[[1]][5]
                                All_data_Heavy[heavy_num,"Total-IMGT_mismatch"]= v_string[[1]][6]
                                All_data_Heavy[heavy_num,"Total-IMGT_gaps"]= v_string[[1]][7]
                                All_data_Heavy[heavy_num,"Total-IMGT_identity_percent"]= v_string[[1]][8]
                                c= c + 1
                            }
                        }
                    }
          
                }
                 #if it is a light chain - kappa and lambda
                if(Chain == "Kappa"|| Chain == "Lambda")
                {
                    light_num = light_num + 1
                    All_data_Light[light_num,"Query"] = gsub("Query= ","",singleton[1])
                    All_data_Light[light_num,"Length"] = gsub("Length=", "", singleton[3])
                    v_string = strsplit(singleton[c], " ")
                    All_data_Light[light_num,"Score_VL"] =v_string[[1]][90]
                    All_data_Light[light_num,"E_value_VL"] = gsub("     ","",v_string[[1]][5])
          
                    while(grepl("IGKV", singleton[c])||grepl("IGHV", singleton[c])||grepl("IGLV", singleton[c])) {c = c + 1}
                    if((grepl("IGKJ", singleton[c]))||(grepl("IGDH", singleton[c])))
                    {
                        v_string = strsplit(singleton[c], " ")
                        All_data_Light[light_num,"Score_JL"] = gsub("</a","",v_string[[1]][4])
                        All_data_Light[light_num,"E_value_JL"] = gsub("     ","",v_string[[1]][5])
                        while((grepl("IGKJ", singleton[c]))||(grepl("IGDH", singleton[c]))) {c = c + 1}
                    }
                    while(!(grepl("Chain\ type", singleton[c]))){c = c + 1 }
                    v_string = strsplit(singleton[c+1], ">")
                    All_data_Light[light_num,"VL_top_match"] = gsub("</td","",v_string[[1]][3])
                    All_data_Light[light_num,"JL_top_match"] = gsub("</td","",v_string[[1]][5])
                    All_data_Light[light_num,"Chain_type"]= gsub("</td","",v_string[[1]][7])
                    All_data_Light[light_num,"Stop_codon"]= gsub("</td","",v_string[[1]][9])
                    All_data_Light[light_num,"V-Jframe"]= gsub("</td","",v_string[[1]][11])
                    All_data_Light[light_num,"Productive"]= gsub("</td","",v_string[[1]][13])
                    All_data_Light[light_num,"Strand"]= gsub("</td","",v_string[[1]][15])
          
                    while(!(grepl("junction\ details", singleton[c]))){c = c + 1 }
                    v_string = strsplit(singleton[c + 3], ">")
                    All_data_Light[light_num,"V_region_end"]= gsub("</td","",v_string[[1]][3])
                    All_data_Light[light_num,"V-J_junction"]= gsub("</td","",v_string[[1]][5])
                    All_data_Light[light_num,"J-region_start"]= gsub("</td","",v_string[[1]][7])
                    while(c != line_count_singleton)
                    {
                        c = c + 1
                        if(grepl("Alignment\ summary",singleton[c]))
                        {
                            c = c + 2
                            if(grepl("FR1-IMGT",singleton[c]))
                            {
                                v_string = strsplit(singleton[c], ">")
                                All_data_Light[light_num,"V_FR1-IMGT_from"]= gsub("</td","",v_string[[1]][5])
                                All_data_Light[light_num,"V_FR1-IMGT_to"]= gsub("</td","",v_string[[1]][7])
                                All_data_Light[light_num,"V_FR1-IMGT_length"]= gsub("</td","",v_string[[1]][9])
                                All_data_Light[light_num,"V_FR1-IMGT_match"]= gsub("</td","",v_string[[1]][11])
                                All_data_Light[light_num,"V_FR1-IMGT_mismatch"]= gsub("</td","",v_string[[1]][13])
                                All_data_Light[light_num,"V_FR1-IMGT_gaps"]= gsub("</td","",v_string[[1]][15])
                                All_data_Light[light_num,"V_FR1-IMGT_identity_percent"]= gsub("</td","",v_string[[1]][17])
                                c= c + 1
                            }
                            if(grepl("CDR1-IMGT",singleton[c]))
                            {
                                v_string = strsplit(singleton[c], ">")
                                All_data_Light[light_num,"CDR1-IMGT_from"]= gsub("</td","",v_string[[1]][5])
                                All_data_Light[light_num,"CDR1-IMGT_to"]= gsub("</td","",v_string[[1]][7])
                                All_data_Light[light_num,"CDR1-IMGT_length"]= gsub("</td","",v_string[[1]][9])
                                All_data_Light[light_num,"CDR1-IMGT_match"]= gsub("</td","",v_string[[1]][11])
                                All_data_Light[light_num,"CDR1-IMGT_mismatch"]= gsub("</td","",v_string[[1]][13])
                                All_data_Light[light_num,"CDR1-IMGT_gaps"]= gsub("</td","",v_string[[1]][15])
                                All_data_Light[light_num,"CDR1-IMGT_identity_percent"]= gsub("</td","",v_string[[1]][17])
                                c= c + 1
                            }
                            if(grepl("FR2-IMGT",singleton[c]))
                            {
                                v_string = strsplit(singleton[c], ">")
                                All_data_Light[light_num,"FR2-IMGT_from"]= gsub("</td","",v_string[[1]][5])
                                All_data_Light[light_num,"FR2-IMGT_to"]= gsub("</td","",v_string[[1]][7])
                                All_data_Light[light_num,"FR2-IMGT_length"]= gsub("</td","",v_string[[1]][9])
                                All_data_Light[light_num,"FR2-IMGT_match"]= gsub("</td","",v_string[[1]][11])
                                All_data_Light[light_num,"FR2-IMGT_mismatch"]= gsub("</td","",v_string[[1]][13])
                                All_data_Light[light_num,"FR2-IMGT_gaps"]= gsub("</td","",v_string[[1]][15])
                                All_data_Light[light_num,"FR2-IMGT_identity_percent"]= gsub("</td","",v_string[[1]][17])
                                c= c + 1
                            }
                            if(grepl("CDR2-IMGT",singleton[c]))
                            {
                                v_string = strsplit(singleton[c], ">")
                                All_data_Light[light_num,"CDR2-IMGT_from"]= gsub("</td","",v_string[[1]][5])
                                All_data_Light[light_num,"CDR2-IMGT_to"]= gsub("</td","",v_string[[1]][7])
                                All_data_Light[light_num,"CDR2-IMGT_length"]= gsub("</td","",v_string[[1]][9])
                                All_data_Light[light_num,"CDR2-IMGT_match"]= gsub("</td","",v_string[[1]][11])
                                All_data_Light[light_num,"CDR2-IMGT_mismatch"]= gsub("</td","",v_string[[1]][13])
                                All_data_Light[light_num,"CDR2-IMGT_gaps"]= gsub("</td","",v_string[[1]][15])
                                All_data_Light[light_num,"CDR2-IMGT_identity_percent"]= gsub("</td","",v_string[[1]][17])
                                c= c + 1
                            }
                            if(grepl("FR3-IMGT",singleton[c]))
                            {
                                v_string = strsplit(singleton[c], ">")
                                All_data_Light[light_num,"FR3-IMGT_from"]= gsub("</td","",v_string[[1]][5])
                                All_data_Light[light_num,"FR3-IMGT_to"]= gsub("</td","",v_string[[1]][7])
                                All_data_Light[light_num,"FR3-IMGT_length"]= gsub("</td","",v_string[[1]][9])
                                All_data_Light[light_num,"FR3-IMGT_match"]= gsub("</td","",v_string[[1]][11])
                                All_data_Light[light_num,"FR3-IMGT_mismatch"]= gsub("</td","",v_string[[1]][13])
                                All_data_Light[light_num,"FR3-IMGT_gaps"]= gsub("</td","",v_string[[1]][15])
                                All_data_Light[light_num,"FR3-IMGT_identity_percent"]= gsub("</td","",v_string[[1]][17])
                                c= c + 1
                                }
                            if(grepl("CDR3-IMGT",singleton[c]))
                            {
                                v_string = strsplit(singleton[c], ">")
                                All_data_Light[light_num,"CDR3-IMGT_germline_from"]= gsub("</td","",v_string[[1]][5])
                                All_data_Light[light_num,"CDR3-IMGT_germline_to"]= gsub("</td","",v_string[[1]][7])
                                All_data_Light[light_num,"CDR3-IMGT_germline_length"]= gsub("</td","",v_string[[1]][9])
                                All_data_Light[light_num,"CDR3-IMGT_germline_match"]= gsub("</td","",v_string[[1]][11])
                                All_data_Light[light_num,"CDR3-IMGT_germline_mismatch"]= gsub("</td","",v_string[[1]][13])
                                All_data_Light[light_num,"CDR3-IMGT_germline_gaps"]= gsub("</td","",v_string[[1]][15])
                                All_data_Light[light_num,"CDR3-IMGT_germline_identity_percent"]= gsub("</td","",v_string[[1]][17])
                                c= c + 1
                            }
                            if(grepl("Total",singleton[c]))
                            {
                                v_string = strsplit(singleton[c], ">")
                                All_data_Light[light_num,"Total-IMGT_length"]= gsub("</td","",v_string[[1]][9])
                                All_data_Light[light_num,"Total-IMGT_match"]= gsub("</td","",v_string[[1]][11])
                                All_data_Light[light_num,"Total-IMGT_mismatch"]= gsub("</td","",v_string[[1]][13])
                                All_data_Light[light_num,"Total-IMGT_gaps"]= gsub("</td","",v_string[[1]][15])
                                All_data_Light[light_num,"Total-IMGT_identity_percent"]= gsub("</td","",v_string[[1]][17])
                                c= c + 1
                            }
                        }
                    }
                }
            }
            rm(c, line_count_singleton, v_string, singleton)
        }
  }
  # save the results table
  All_data_Heavy <- All_data_Heavy[rowSums(is.na(All_data_Heavy)) != ncol(All_data_Heavy),]
  if(length(All_data_Heavy) != 0)
  {
    file_heavy = paste("./Parse_IgBlast/","Parse","_heavy_", gsub("Igblast_results_","",file_name), sep = "")
    write.table(All_data_Heavy, file_heavy, sep="\t")
  }
  
  All_data_Light <- All_data_Light[rowSums(is.na(All_data_Light)) != ncol(All_data_Light),]
  if(length(All_data_Light) != 0)
  {
    file_light = paste("./Parse_IgBlast/","Parse","_light_", gsub("Igblast_results_","",file_name), sep = "")
    write.table(All_data_Light, file_light, sep="\t")
  }
  rm(All_data_Light,All_data_Heavy, file_name, Chain,condition,heavy_num,light_num, file_heavy, long,main_file_line_count,seq_number,number_of_sequences,rawHTML)
}


setwd("C:/Users/Chen-Lab/Desktop/Program/BCR_rep/IgBlast_result")
# find all the txt file
seqs <- Sys.glob("*.txt")
# loop all the file with function ParseIgBlast
for (i in 1:length(seqs)){
  start <- Sys.time()
  print(seqs[i])
  ParseIgBlast(seqs[i])
  end <- Sys.time()
  print(paste0("运行时间是：", (end-start)))
}

rm(i, seq, seqs, end, start)