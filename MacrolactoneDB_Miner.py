#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      zinph
#
# Created:     01/11/2016
# Copyright:   (c) zinph 2016
# Licence:     <your licence>
#-------------------------------------------------------------------------------


import pandas as pd
from rdkit import RDConfig
from rdkit.Chem import PandasTools
import time
import pickle

'''
RS_min = int
RS_max = int
MW_min = int
MW_max = int
nRing_min = int
nRing_max = int
Lipinski = 'y' or 'n'
nG12Ring_min = int
nG12Ring_max = int
SlogP_min
SlogP_max
nSugars_min
nSugars_max
core_ester_min
core_ester_max
naRing_min  # number of aromatic ring count
naRing_max
'''

class MacrolactoneDB_Miner:

    def __init__(self,command):
        self.df = pd.read_csv('MacrolactoneDB.csv')
        self.command = command

    def limit_RS(self, df, min_RS, max_RS):
        '''
        Limit the pandas frame by ring size and return the filtered inchi keys.
        '''
        if min_RS=='dc' and max_RS=='dc':
            df = df
        else:
            if min_RS!='dc' and max_RS!='dc':
                df = df[(df['SmallestRS'] >= min_RS) & (df['LargestRS'] <= max_RS)]
            elif min_RS=='dc' and max_RS:
                df = df[df['LargestRS'] <= max_RS]
            elif max_RS=='dc' and min_RS:
                df = df[df['SmallestRS'] >= min_RS]
        return set(df['InChI Keys'].unique())

    def limit_SlogP(self, df, min_SlogP, max_SlogP):
        '''
        Limit the pandas frame by ring size and return the filtered inchi keys.
        '''
        if min_SlogP=='dc' and max_SlogP=='dc':
            df = df
        else:
            if min_SlogP!='dc' and max_SlogP!='dc':
                df = df[(df['SLogP'] >= min_SlogP) & (df['SLogP'] <= max_SlogP)]
            elif min_SlogP=='dc' and max_SlogP:
                df = df[df['SLogP'] <= max_SlogP]
            elif max_SlogP=='dc' and min_SlogP:
                df = df[df['SLogP'] >= min_SlogP]
        return set(df['InChI Keys'].unique())

    def limit_core_ester(self, df, min_core_ester, max_core_ester):
        '''
        Limit the pandas frame by ring size and return the filtered inchi keys.
        '''
        if min_core_ester=='dc' and max_core_ester=='dc':
            df = df
        else:
            if min_core_ester!='dc' and max_core_ester!='dc':
                df = df[(df['core_ester'] >= min_core_ester) & (df['core_ester'] <= max_core_ester)]
            elif min_core_ester=='dc' and max_core_ester:
                df = df[df['core_ester'] <= max_core_ester]
            elif max_core_ester=='dc' and min_core_ester:
                df = df[df['core_ester'] >= min_core_ester]
        return set(df['InChI Keys'].unique())

    def limit_MW(self, df, min_MW, max_MW):
        '''
        Limit the pandas frame by molecular weight and return the filtered inchi keys.
        '''
        if min_MW=='dc' and max_MW=='dc':
            df = df
        else:
            if min_MW!='dc' and max_MW!='dc':
                df = df[(df['MW'] >= min_MW) & (df['MW'] <= max_MW)]
            elif min_MW=='dc' and max_MW:
                df = df[df['MW'] <= max_MW]
            elif max_MW=='dc' and min_MW:
                df = df[df['MW'] >= min_MW]
        return set(df['InChI Keys'].unique())

    def limit_nRing(self, df, min_nRing, max_nRing):
        '''
        Limit the pandas frame by number of rings and return the filtered inchi keys.
        '''
        if min_nRing=='dc' and max_nRing=='dc':
            df = df
        else:
            if min_nRing!='dc' and max_nRing!='dc':
                df = df[(df['nRing'] >= min_nRing) & (df['nRing'] <= max_nRing)]
            elif min_nRing=='dc' and max_nRing:
                df = df[df['nRing'] <= max_nRing]
            elif max_nRing=='dc' and min_nRing:
                df = df[df['nRing'] >= min_nRing]
        return set(df['InChI Keys'].unique())

    def limit_Lipinski(self, df, Ro5='n'):
        Ro5 = Ro5.lower()
        if Ro5=='dc':
            df = df
        else:
            if Ro5 == 'TRUE' or Ro5 == 'y' or Ro5 == 'yes':
                df = df[df['Lipinski'].astype(str).str.contains('T')]
            elif Ro5 == 'FALSE' or Ro5 == 'n' or Ro5 == 'no':
                df = df[df['Lipinski'].astype(str).str.contains('F')]
        return set(df['InChI Keys'].unique())

    def limit_nFusedRing(self, df, min_nFRing, max_nFRing):
        if min_nFRing=='dc' and max_nFRing=='dc':
            df = df
        else:
            if min_nFRing!='dc' and max_nFRing!='dc':
                df = df[(df['nFRing'] >= min_nFRing) & (df['nFRing'] <= max_nFRing)]
            elif min_nFRing=='dc' and max_nFRing:
                df = df[df['nFRing'] <= max_nFRing]
            elif max_nFRing=='dc' and min_nFRing:
                df = df[df['nFRing'] >= min_nFRing]
        return set(df['InChI Keys'].unique())


    def limit_nG12Ring(self, df, min_nG12Ring, max_nG12Ring):
        if min_nG12Ring=='dc' and max_nG12Ring=='dc':
            df = df
        else:
            if min_nG12Ring!='dc' and max_nG12Ring!='dc':
                df = df[(df['nG12Ring'] >= min_nG12Ring) & (df['nG12Ring'] <= max_nG12Ring)]
            elif min_nG12Ring=='dc' and max_nG12Ring:
                df = df[df['nG12Ring'] <= max_nG12Ring]
            elif max_nG12Ring=='dc' and min_nG12Ring:
                df = df[df['nG12Ring'] >= min_nG12Ring]
        return set(df['InChI Keys'].unique())


    def limit_nSugars(self, df, min_nSugars, max_nSugars):
        if min_nSugars=='dc' and max_nSugars=='dc':
            df = df
        else:
            if min_nSugars!='dc' and max_nSugars!='dc':
                df = df[(df['nSugars'] >= min_nSugars) & (df['nSugars'] <= max_nSugars)]
            elif min_nSugars=='dc' and max_nSugars:
                df = df[df['nSugars'] <= max_nSugars]
            elif max_nSugars=='dc' and min_nSugars:
                df = df[df['nSugars'] >= min_nSugars]
        return set(df['InChI Keys'].unique())

    def limit_naRing(self, df, min_naRing, max_naRing):
        '''
        Limit the pandas frame by ring size and return the filtered inchi keys.
        '''
        if min_naRing=='dc' and max_naRing=='dc':
            df = df
        else:
            if min_naRing!='dc' and max_naRing!='dc':
                df = df[(df['naRing'] >= min_naRing) & (df['naRing'] <= max_naRing)]
            elif min_naRing=='dc' and max_naRing:
                df = df[df['naRing'] <= max_naRing]
            elif max_naRing=='dc' and min_naRing:
                df = df[df['naRing'] >= min_naRing]
        return set(df['InChI Keys'].unique())

    def limit_activity_reported(self, df, activity='dc'):
        if activity.lower() == 'yes' or activity.lower() == 'y' or activity == 1:
            df = df[df['activity_reported'].astype(str).str.contains('yes')]
        elif activity.lower() =='no' or activity.lower() =='n' or activity == 0:
            df = df[df['activity_reported'].astype(str).str.contains('yes')]
        else:
            df = df
        return set(df['InChI Keys'].unique())

    def frame_manage(self):
        PandasTools.AddMoleculeColumnToFrame(self.filtered_df,'smiles','structures')
        structures = self.filtered_df['structures']
        df = self.filtered_df.drop(columns=['structures'])
        df.insert(0, 'structures', structures)
        # df = df[['ChEMBL_IDs','structures',"target_organism","target_molecule_pref_name",'# Known Targets','Known Targets','smiles']]
        df = df[['IDs', 'structures', 'molecule_pref_name', '# Target Organisms', 'Target Organisms', '# Known Targets', 'Known Targets', 'target_pref_name']]
        df = df.sort_values(by='# Known Targets', ascending=False, na_position='last')
        return df

    def compile_filters(self):
        RS_inchi = self.limit_RS(self.df, self.command['RS_min'], self.command['RS_max'])
        MW_inchi = self.limit_MW(self.df, self.command['MW_min'], self.command['MW_max'])
        nRing_inchi = self.limit_nRing(self.df, self.command['nRing_min'], self.command['nRing_max'])
        Lipinski_inchi = self.limit_Lipinski(self.df, self.command['Lipinski'])
        nG12Ring_inchi = self.limit_nG12Ring(self.df, self.command['nG12Ring_min'], self.command['nG12Ring_max'])
        SlogP_inchi = self.limit_SlogP(self.df, self.command['SlogP_min'], self.command['SlogP_max'])
        Sugars_inchi = self.limit_nSugars(self.df, self.command['nSugars_min'], self.command['nSugars_min'])
        nFRing_inchi = self.limit_nFusedRing(self.df, self.command['nFRing_min'], self.command['nFRing_max'])
        core_ester_inchi = self.limit_core_ester(self.df, self.command['core_ester_min'], self.command['core_ester_max'])
        naRing_inchi = self.limit_naRing(self.df, self.command['naRing_min'], self.command['naRing_max'])
        activity_reported_inchi = self.limit_activity_reported(self.df, self.command['activity_reported'])

        sets = [RS_inchi, MW_inchi, nRing_inchi, Lipinski_inchi, nG12Ring_inchi,SlogP_inchi,Sugars_inchi,nFRing_inchi,core_ester_inchi,naRing_inchi,activity_reported_inchi]
        self.filtered_inchi = list(set.intersection(*sets))
        self.filtered_df = self.df.loc[self.df['InChI Keys'].isin(self.filtered_inchi)]
        # print(filtered_df.shape[0], ' compouds have been compiled based on your filters.')
        # smiles = filtered_df['smiles'].tolist()
        PandasTools.AddMoleculeColumnToFrame(self.filtered_df,'smiles','Molecule picture')

        # export csv file
        # self.filtered_df.to_csv('temp.csv', index=False)

        ## export sdf file
        # PandasTools.WriteSDF(self.filtered_df, 'temp.sdf', molColName='structures', properties=list(self.filtered_df.columns), allNumeric=False)

        # export smiles
        # self.smiles_writer()

        # self.filtered_df.to_sql(name='temp', con=db.engine, index=False)
        smiles_frame = self.frame_manage()

        return smiles_frame
        # filtered_df['smiles'].to_csv(filename,index=False)

    def inchi_writer(self):
        f = open('filtered_inchi','w')
        f.write(','.join(self.filtered_df['InChI Keys'].tolist()))
        f.close()

    def smiles_writer(self):
        f = open('temp.smiles','w')
        f.write('smiles\n')
        f.write('\n'.join(self.filtered_df['smiles'].tolist()))
        f.close()


def convert_time(second):
    day = second/86400
    hour = (day - int(day))*24
    minute = (hour - int(hour))*60
    second = round((minute - int(minute))*60,4)
    return(str(int(day)) + ' DAYS: '+ str(int(hour)) + ' HOURS: '+ str(int(minute)) + ' MINUTES: ' + str(second) + ' SECONDS')


# def main():
#     command = {
#             'RS_min' : 12,
#             'RS_max' : 16,
#             'MW_min' : 'dc',
#             'MW_max' : 'dc',
#             'nRing_min' : 'dc',
#             'nRing_max' : 'dc',
#             'Lipinski' : 'dc',
#             'nG12Ring_min': 1,
#             'nG12Ring_max' : 1,
#             'SlogP_min' : 'dc',
#             'SlogP_max' : 'dc',
#             'nSugars_min' : 2,
#             'nSugars_max' : 2,
#             'nFRing_min' : 0,
#             'nFRing_max' : 0,
#             'core_ester_min':1,
#             'core_ester_max':1,
#             'naRing_min':0,
#             'naRing_max':0
#             }
#
#     start_time = time.time()
#     sample = MacrolactoneDB(command)
#     sample.compile_filters('classic_macrolides.smiles')
#     duration = convert_time(time.time()-start_time)
#     print('Time Elapsed:' , duration)
# main()
