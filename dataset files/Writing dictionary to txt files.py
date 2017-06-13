import os
central_directory = "C:\\Users\\sanja\\University of Toronto\\OneDrive - University of Toronto\\University of Toronto Medical School\\Year 2 summer\\Summer project\\MitoTIP\\Dictionary files"
os.chdir(central_directory)
import ast
file = open ("dict_complement_pos.txt", "w")
for key in list(dict_complement_pos.keys()): 
    file.write(str(key) +"\n")
    file.write(str(dict_complement_pos[key]) +"\n")
file.close()

file = open ("dict_tRNA_features_stored.txt", "w")
for key in list(dict_tRNA_features_stored.keys()): 
    file.write(str(key) +"\n")
    file.write(str(dict_tRNA_features_stored[key]) +"\n")
file.close()

file = open("dict_complement_pos.txt", "r")
dict_complement_pos_new = {}
for line in file: 
    dict_complement_pos_new[line.strip("\n")] = ast.literal_eval(file.readline().strip("\n"))
    
file = open("dict_tRNA_features_stored.txt", "r")
dict_tRNA_features_stored_new = {}
for line in file: 
    dict_tRNA_features_stored_new[line.strip("\n")] = ast.literal_eval(file.readline().strip("\n"))
    



# dict_complement_pos_stored = dict_complement_pos


# dict_tRNA_features_stored = dict_tRNA_features
