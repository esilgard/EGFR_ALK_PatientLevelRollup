# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 16:23:45 2016

@author: esilgard
"""
import os
## patient/case level roll up rules for egfr/alk classification
## directory structure is based on the output from the report level output
## from the original EGFR_ALK_Classification

run = 'System'
gold_case_labels = {}
report_d = {}

def get_system_report_labels():
    for batch in ['InternalValidationOutput']:
        result_directory = 'H:/EGFR_ALK/' + batch + '/Result/'
        method_directory = 'H:/EGFR_ALK/' + batch + '/Method/'
        method_error_directory = 'H:/EGFR_ALK/' + batch + '/MethodErrors/'
        result_error_directory = 'H:/EGFR_ALK/' + batch + '/ResultErrors/'
        insuff_unk_directory = 'H:/EGFR_ALK/' + batch + '/InsufficientOrUnknown/'
        algorithm_directory_mapping = {result_directory:'Result', method_directory:'Method', insuff_unk_directory:'Result'}
        
        algorithm_directory_mapping[result_error_directory] = 'Result'
        algorithm_directory_mapping[method_error_directory] = 'Method'
    
        for folder,algorithm_string in algorithm_directory_mapping.items():
            for path, directory, files in os.walk(folder):
                for f in files:                
                    instances = [x.strip() for x in open(path+'/'+f,'r').readlines()]
                    headers = dict((k, v) for v, k in enumerate(instances[0].split('\t')))
                    # CaseId      ReportId        TestInstance      InstanceId      SystemOutput  GoldStandardLabel
                                   
                    for inst in instances[1:]:
                        i = inst.split('\t')
                        case = i[headers.get('CaseId')]
                        instance_id = i[headers.get('InstanceId')]
                        test = i[headers.get('TestInstance')]
                        report_d[case] = report_d.get(case,{})
                        report_d[case]['System'] = report_d[case].get('System',{})
                        report_d[case]['System'][test] = report_d[case]['System'].get(test,{})
                        report_d[case]['System'][test][instance_id] = report_d[case]['System'][test].get(instance_id,{})
                        report_d[case]['System'][test][instance_id][algorithm_string] = i[headers.get('SystemOutput')]
                                                                    
def get_gold_report_labels():
    gold_annotations = [x.strip().split('\t') for x in open('H:/EGFR_ALK/AnnotationConsolidation/reformatted_gold_case_inst_labels.txt','r').readlines()[1:]]
    for lines in gold_annotations:
        # CaseId       EGFR CaseLabel     ALK CaseLabel   ReportId   EGFR TestDone   ALK TestDone    EGFR Result  ALK Result  EGFR Method ALK Method
        case = lines[0]                
        if case not in gold_case_labels: 
            gold_case_labels[case] = {}
            gold_case_labels[case]['EGFR'] = tuple(lines[1].split(';'))
            gold_case_labels[case]['ALK'] = tuple(lines[2].split(';'))
        inst = lines[3]
        report_d[case] = report_d.get(case,{})
        report_d[case]['Gold'] = report_d[case].get('Gold',{})
        report_d[case]['Gold']['EGFR'] = report_d[case]['Gold'].get('EGFR',{})
        report_d[case]['Gold']['ALK'] = report_d[case]['Gold'].get('ALK',{})
        report_d[case]['Gold']['EGFR'][inst] = report_d[case]['Gold']['EGFR'].get(inst,{})
        report_d[case]['Gold']['ALK'][inst] = report_d[case]['Gold']['ALK'].get(inst,{})
        
        egfr_result = lines[4]
        if 'Results Reported' in egfr_result:
            egfr_result = lines[6]
        egfr_method = lines[8]
      
        alk_result = lines[5]
        if 'Results Reported' in alk_result:
            alk_result = lines[7]
        alk_method = lines[9]
        
        report_d[case]['Gold']['EGFR'][inst]['Result'] = egfr_result
        report_d[case]['Gold']['EGFR'][inst]['Method'] = egfr_method
        report_d[case]['Gold']['ALK'][inst]['Result'] = alk_result
        report_d[case]['Gold']['ALK'][inst]['Method'] = alk_method



def get_roll_up():
    match = 0
    no_match = 0
    total = 0.0
    case_label = {}
    with open('H:/EGFR_ALK/InternalValidationOutput/'+ run + '_case_level_labels.txt','w') as out:
        out.write('CaseId\tTest\tSystem Label\tGold Label\n')
        for case in report_d:
            case_label[case] = {}
            trumping_roll_up_rules_sort = {'EGFR':[('Positive','MutationalAnalysis'),('Negative','MutationalAnalysis'),('Positive','OTHER'),\
                                                ('Negative','OTHER'),('Insufficient','None'),('Unknown','None'),('NotReported','None')],\
                                            'ALK':[('Positive','FISH'),('Negative','FISH'),('Positive','OTHER'),\
                                                ('Negative','OTHER'),('Insufficient','None'),('Unknown','None'),('NotReported','None')]}
            
            for test_name in ['EGFR','ALK']:
    
                test_report_d = report_d[case][run].get(test_name)
    
                ## default to Unknown/NotReported for the case label
                result_list = [(test_report_d[inst_id].get('Result','Unknown'), test_report_d[inst_id].get('Method','None')) for inst_id in test_report_d]
                # for error propogation mismatches (there's a method, but no result, etc) - assume the most common
                # based on the best performing algorithm (result - pos or neg)
       
                for r in range(len(result_list)):
                    system_result = result_list[r][0]
               
                    if result_list[r] not in trumping_roll_up_rules_sort[test_name]:
                        if system_result == 'Positive': result_list[r] = trumping_roll_up_rules_sort[test_name][0]
                        elif  system_result == 'Negative': result_list[r] = trumping_roll_up_rules_sort[test_name][1]
                        else:
                            result_list[r] = (system_result,'None')
                  
                ## arg max (min really, since the sort is counted from zero ascending) of result sets based on trumping_roll_up_rules_sort 
                best_sorted_result = min([trumping_roll_up_rules_sort[test_name].index(result_tuple) for result_tuple in result_list]) 
    
                case_label[case][test_name] = trumping_roll_up_rules_sort[test_name][best_sorted_result]
                
                if case_label[case][test_name] == gold_case_labels[case][test_name]: match +=1
                else: 
                    no_match +=1    
                out.write(case + '\t' + test_name + '\t' + ' by '.join(case_label[case][test_name]) +  '\t' + ' by '.join(gold_case_labels[case][test_name]) + '\n')
                total +=1
            
    print (str(match/total) + ' matched')
    print (str(no_match/total) + ' error rate')
                

    
                     
get_system_report_labels()
get_gold_report_labels()

print (report_d['PAT-00705030_2'].keys())
patient_labels = get_roll_up()


    
