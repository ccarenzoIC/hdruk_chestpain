from pandas import DataFrame, Series, read_csv
from json import loads
from urllib.request import urlopen
from flatten_json import flatten 
from black import format_file_contents, FileMode
import sys
from requests import get
from bs4 import BeautifulSoup

def generate_dd_pheno(mappings_icd10_snomed_path):

    """
    This function is used to obtain the codes from the phenotypes listed at the protocol producing a dictionary
    with the category, the subcategory and the mapping systems used.
   
    """
    phenotypes_dictionary = {
			"hypertension": "http://phenotypes.healthdatagateway.org/api/v1/public/phenotypes/PH189/version/378/export/codes/?format=json",
			"smoking": "http://phenotypes.healthdatagateway.org/api/v1/public/phenotypes/PH982/version/2160/export/codes/?format=json",
			"diabetes": "http://phenotypes.healthdatagateway.org/api/v1/public/phenotypes/PH152/version/304/export/codes/?format=json",
			"myocardial_infarction": "http://phenotypes.healthdatagateway.org/api/v1/public/phenotypes/PH215/version/430/export/codes/?format=json",
			"heart_failure": "http://phenotypes.healthdatagateway.org/api/v1/public/phenotypes/PH182/version/364/export/codes/?format=json",
			"stroke": "http://phenotypes.healthdatagateway.org/api/v1/public/phenotypes/PH85/version/170/export/codes/?format=json",
			"ischaemic_stroke": "http://phenotypes.healthdatagateway.org/api/v1/public/phenotypes/PH56/version/112/export/codes/?format=json"
			}

    #TODO: Use API for mapping instead of static file.
    mapping_icd10_snomed = read_csv(mappings_icd10_snomed_path)
    out = {}
    aux_result = []
    for k, v in phenotypes_dictionary.items():
            response = urlopen(v)
            data_json = loads(response.read())
            df = DataFrame(data_json)
            aux_dic_list = []
            aux_dic = {}
            if(df[df.coding_system == "ICD10 codes"].size):
                #getting ICD10 codes, df preparation
                df = df[df.coding_system == "ICD10 codes"]
                df["code"] = df.code.str.replace(".","",regex = True)
                df = df.rename(columns = {"code":"icd10_code","concept_id":"concept_id_ph","concept_version_id":"concept_version_id_ph","description":"description_ph"})
                df[["disease","category"]]=df['code_attributes'].apply(flatten).apply(Series)
                # merging icd10, snomed df with icd10, phenotype
                mapping_icd10_snomed["icd10_term"] = mapping_icd10_snomed.icd10_code.apply(lambda x: x[0:3])
                df_phenotype = mapping_icd10_snomed.merge(df[["icd10_code","category","concept_id_ph","concept_version_id_ph","phenotype_id","phenotype_version_id","phenotype_name","disease", "description_ph"]], on = "icd10_code")    
                if(df_phenotype.size == 0):
                    df_phenotype = mapping_icd10_snomed.merge(df[["icd10_code","category","concept_id_ph","concept_version_id_ph","phenotype_id","phenotype_version_id","phenotype_name","disease", "description_ph"]], left_on = "icd10_term" ,right_on = "icd10_code")         
                    df_phenotype = df_phenotype.rename(columns = {"icd10_term":"icd10_code"})
                for subcat in df_phenotype.category.drop_duplicates().values:
                    members_icd10 = [
                    i
                    for i in df_phenotype[df_phenotype.category==subcat].icd10_code.astype(str).drop_duplicates().values
                    ]
                    members_snomed = [
                    int(i)
                    for i in df_phenotype[df_phenotype.category==subcat].snomed_concept_id.astype(str).astype(float).drop_duplicates().values
                    ]
                    aux_result = [{"mapping_system":"icd_10","members":members_icd10},{"mapping_system":"snomed","members":members_snomed}]
                    aux_dic[subcat] = aux_result
                    aux_dic_list.append(aux_dic)
                    aux_dic = {}
                    aux_result = []
            elif(df[df.coding_system == "SNOMED  CT codes"].size):
                df = df.rename(columns = {"code":"snomed_concept_id","concept_id":"concept_id_ph","concept_version_id":"concept_version_id_ph","description":"description_ph"})
                df["disease"] = ""
                df["category"] = df["concept_name"]
                mapping_icd10_snomed.snomed_concept_id = mapping_icd10_snomed.snomed_concept_id.astype(str)
                df_phenotype = mapping_icd10_snomed.merge(df[["snomed_concept_id","category","concept_id_ph","concept_version_id_ph","phenotype_id","phenotype_version_id","phenotype_name","disease", "description_ph"]], on = "snomed_concept_id")
                for subcat in df_phenotype.category.drop_duplicates().values:
                    members_icd10 = [
                    i
                    for i in df_phenotype[df_phenotype.category==subcat].icd10_code.astype(str).drop_duplicates().values
                    ]
                    members_snomed = [
                    int(i)
                    for i in df_phenotype[df_phenotype.category==subcat].snomed_concept_id.astype(str).astype(float).drop_duplicates().values
                    ]
                    aux_result = [{"mapping_system":"icd_10","members":members_icd10},{"mapping_system":"snomed","members":members_snomed}]
                    aux_dic[subcat] = aux_result
                    aux_dic_list.append(aux_dic)
                    aux_dic = {}
                    aux_result = []
            out[k] = {"url": v, "mapping": aux_dic_list}    
    
    dd_py = ""
    for k, v in out.items():
        line = f"{k} = {v}\n\n"
        # replace non-ascii characters
        line = line.encode(encoding="ascii", errors="replace").decode().replace("?", "")
        #line = line.replace('(','').replace(')','')  
        dd_py += line        

    dd_py = format_file_contents(dd_py, fast=True, mode=FileMode())
    with open("./chest_pain/data/dd.py", "wt") as f:
        f.write(dd_py)

    print("Successfully created black formatted dd.py")

    print(dd_py)
        
    
def generate_dd_opcs(mappings_snomed_opcs_path):
    
    """
    This function is used to obtain the mapping to snomed from the OPCS codes listed at the protocol producting a dictionary
    with the concept, the codes and the mapping systems used.
   
    """

    opcs_dictionary = {
			"pci_procedure" : ["K49", "K49.1", "K49.2", "K49.3","K49.4", "K49.8", "K49.9", "K50",
								"K50.1", "K50.4", "K50.8", "K50.9","K75", "K75.1", "K75.2", "K75.3", 
								"K75.4", "K75.8", "K75.9"],
			"coronary_bypass_grafting" : ["K49", "K49.1", "K49.2", "K49.3","K49.4", "K49.8", "K49.9", "K50",
											"K50.1", "K50.4", "K50.8", "K50.9","K75", "K75.1", "K75.2", "K75.3", 
											"K75.4", "K75.8", "K75.9"],
            "angiography" : ["K63.1", "K63.2", "K63.3", "K63.4", "K63.5", "K63.6", "K65.1", "K65.2", "K65.3", "K65.8", "K65.9" ]
			}
    #TODO: Use API for mapping instead of static file.
    mapping_opcs_snomed = read_csv(mappings_snomed_opcs_path)
    out = {}
    aux_dic_list = []
    for k, v in opcs_dictionary.items():
        aux_dic = {}
        aux = DataFrame(v, columns = ["opcs_term"])
        aux.opcs_term = aux.opcs_term.str.replace(".","",regex=True)
        aux_merge = aux.merge(mapping_opcs_snomed[["opcs_term","snomed_concept_id"]], on = "opcs_term")
        members_opcs = [
        i
        for i in aux_merge.opcs_term.drop_duplicates().values
        ]
        members_snomed = [
        int(i)
        for i in aux_merge.snomed_concept_id.drop_duplicates().values
        ]
        aux_dic = [{"mapping_system":"icd_10","members":members_opcs},{"mapping_system":"snomed","members":members_snomed}]
        aux_dic_list.append(aux_dic)

        out[k] = {"mapping": aux_dic_list}

    dd_opcs = ""
    for k, v in out.items():
        line = f"{k} = {v}\n\n"
        # replace non-ascii characters
        line = line.encode(encoding="ascii", errors="replace").decode().replace("?", "")
        dd_opcs += line        

    dd_opcs = format_file_contents(dd_opcs, fast=True, mode=FileMode())
    with open("./chest_pain/data/dd_opcs.py", "wt") as f:
        f.write(dd_opcs)

    print("Successfully created black formatted dd_opcs.py")

    print(dd_opcs)

def generate_dd_icd10(mappings_icd10_snomed_path):
    """
    This function is used to obtain the mapping to snomed from the ICD10 codes listed at the protocol producting a dictionary
    with the concept, the codes and the mapping systems used.    
    """
    icd10_dictionary = {
			  "cardiovascular_death" : ["I00","I01","I02","I03","I04","I05","I06","I07","I08","I09","I10","I11","I12",
                             		"I13","I14","I15","I16","I17","I18","I19","I20","I21","I22","I23","I24","I25","I26","I27",
									"I28","I29","I30","I31","I32","I33","I34","I35","I36","I37","I38","I39","I40","I41","I42","I43",
									"I44","I45","I46","I47","I48","I49","I50","I51","I52","I53","I54","I55","I56","I57","I58","I59",
									"I60","I61","I62","I63","I64","I65","I66","I67","I68","I69","I70","I71","I72","I73","I74","I75",
                                    "I76","I77","I78","I79","I80","I81","I82","I83","I84","I85","I86","I87","I88","I89","I90","I91",
									"I92","I93","I94","I95","I96","I97","I98","I99"],
                "cardiac_death" : ["I05","I06","I07","I08","I09","I20","I21","I22","I23","I24","I25","I30","I31","I32","I33","I34",
                                   "I35","I36","I37","I38","I39","I40","I41","I42","I43","I44","I45","I46","I47","I48","I49","I50","I51"]
				}

    #TODO: Use API for mapping instead of static file.
    mapping_icd10_snomed = read_csv(mappings_icd10_snomed_path)
    out = {}
    aux_dic = {}
    aux_dic_list = []
    for k, v in icd10_dictionary.items():
        aux = DataFrame(v, columns = ["icd10_code"])
        aux_merge = aux.merge(mapping_icd10_snomed[["icd10_code","snomed_concept_id"]], on = "icd10_code")
        mapping_icd10_snomed["icd10_term"] = mapping_icd10_snomed.icd10_code.apply(lambda x: x[0:3])
        if(aux_merge.size == 0):
            aux_merge = aux.merge(mapping_icd10_snomed[["icd10_code","snomed_concept_id", "icd10_term"]],left_on = "icd10_code", right_on = "icd10_term" )         
            aux_merge = aux_merge.rename(columns = {"icd10_term":"icd10_code"})
        members_icd10 = [
        i
        for i in aux_merge.icd10_code.drop_duplicates().values
        ]
        members_snomed = [
        int(i)
        for i in aux_merge.snomed_concept_id.drop_duplicates().values
        ]
        aux_dic = [{"mapping_system":"icd_10","members":members_icd10},{"mapping_system":"snomed","members":members_snomed}]
        aux_dic_list.append(aux_dic)
        out[k] = {"mapping": aux_dic_list}

    dd_icd10 = ""
    for k, v in out.items():
        line = f"{k} = {v}\n\n"
        # replace non-ascii characters
        line = line.encode(encoding="ascii", errors="replace").decode().replace("?", "")
        dd_icd10 += line        

    dd_icd10 = format_file_contents(dd_icd10, fast=True, mode=FileMode())
    with open("./chest_pain/data/dd_icd10.py", "wt") as f:
        f.write(dd_icd10)

    print("Successfully created black formatted dd_icd10.py")

    print(dd_icd10)


    def generate_dd_bnf():
  
        bnf_dictionary = {
                        "lipid_regulating_drugs":{"description":"BNF Chapter 2.12","url":"https://openprescribing.net/bnf/0212/"},
                        "anti_platelet_drugs":{"description":"BNF chapter 2.9","url":"https://openprescribing.net/bnf/0209/"},
                        "beta_adrenoceptor_blocking_drugs":{"description":"BNF Chapter 2.4","url":"https://openprescribing.net/bnf/0204/"},
                        "renin_angiotensin_system_drugs":{"description":"BNF chapter 2.5.5","url":"https://openprescribing.net/bnf/020505/"}
                        }


        # Get SNOMED Codes and member codes
        drugs_members = {}
        for k, v in bnf_dictionary.items():
            print("Getting list of drugs for ", k)
            r = get(v["url"])
            soup = BeautifulSoup(r.text, 'html.parser')
            mainview = soup.find(class_="starter-template")
            out = []
            for child in mainview.children:
                if child.name == 'a':
                    drug = child.text
                    drug = drug[0:drug.find("(")]
                    out.append(drug)

            print("Got", len(out), " drugs")
            drugs_members[k] = {"description":v["description"],"drug_list":out}

        # Convert dict object into a python module
        text = ""
        for k, v in drugs_members.items():
            line = f"{k} = {v}\n\n"
            # replace non-ascii characters
            line = line.encode(encoding="ascii", errors="replace").decode().replace("?", "")
            text += line
        print(text)
        text = format_file_contents(text, fast=True, mode=FileMode())

        with open("./chest_pain/data/dd_bnf.py", "wt") as f:
            f.write(text)

    if __name__ == "__main__":
        # generate_nhsdd()
        # generate_nhsdd_snomed()
        if sys.argv[1] == "pheno":
            print("Generating chest_pain/data/dd.py")
            generate_dd_pheno("./data/interim/mappings_icd10.csv")
        elif sys.argv[1] == "icd10":
            print("Generating chest_pain/data/dd_icd10.py")
            generate_dd_icd10("./data/interim/mappings_icd10.csv")
        elif sys.argv[1] == "opcs":
            print("Generating avoidable_admissions/data/dd_opcs.py")
            generate_dd_opcs("./data/interim/mappings_opcs.csv")
        else:
            print(
                "Generating chest_pain/data/dd.py, chest_pain/data/dd_icd10.py and chest_pain/data/dd_opcs.py"
            )
            generate_dd_pheno(".data/interim/mappings_icd10.csv")
            generate_dd_icd10(".data/interim/mappings_icd10.csv")
            generate_dd_opcs(".data/interim/mappings_opcs.csv")

            