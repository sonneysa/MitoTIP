# for pos_num in sorted(list(mut_dict.keys())):
#     if mut_dict[pos_num][0]>1 and mut_dict[pos_num][1]>1:
#         print (str(pos_num) + ": " + str(mut_dict[pos_num]))
''' Only updated column I of first sheet'''


central_directory = "C:\\Users\\sanja\\University of Toronto\\OneDrive - University of Toronto\\University of Toronto Medical School\\Year 2 summer\\Summer project\\MitoTIP\\Final MitoTIP uptodate"
dataset = "\\dataset files"
output_dir = "\\Output"
os.chdir(central_directory + dataset)
workbook_name = "PON-mt dataset"
wb = load_workbook(workbook_name + '.xlsx', data_only = True)
PON_dataset = wb['PON database']
list_PON_mut = []
list_PON_poly = []
dict_to_recalc = {'Polymorphic': 'poly', 'Pathologic': 'path'}
output_PON_predictive = {}
output_PON_predictive['path'] = {}
output_PON_predictive['poly'] = {}
classified_mut = 0
classified_poly = 0
threshold = 14.3

for i in range (5,151): #5 to 150 inclusive 
    pos_num = PON_dataset.cell(row = i, column = 1).value 
    base = PON_dataset.cell(row = i, column = 3).value 
    # recalc = dict_to_recalc[PON_dataset.cell(row = i, column = 8).value]
    recalc = None
    output = get_predictive_score(pos_num, base, recalc)
    recalc = dict_to_recalc[PON_dataset.cell(row = i, column = 8).value]
    if recalc == 'path':
        list_PON_mut.append([pos_num, base])
        if output[1][3]>threshold:
            classified_mut += 1
        if output_PON_predictive['path'].get(pos_num, None) == None:
            output_PON_predictive['path'][pos_num] = {}
        output_PON_predictive['path'][pos_num][base] = output[1][3]
    elif recalc == 'poly':
        list_PON_poly.append([pos_num, base])
        if output[1][3]<=threshold:
            classified_poly += 1
        if output_PON_predictive['path'].get(pos_num, None) == None:
            output_PON_predictive['path'][pos_num] = {}
        output_PON_predictive['path'][pos_num][base] = output[1][3]
    PON_dataset.cell(row=i, column = 9, value = output[1][3])
    # PON_dataset.cell(row=i, column = 9, value = 5)
os.chdir(central_directory + output_dir)
wb.save("PON-mt dataset_scored" + '.xlsx')
print ("Sensitivity:" + str(classified_mut/len(list_PON_mut)))
print("Specificity: " + str(classified_poly/len(list_PON_poly)))
    