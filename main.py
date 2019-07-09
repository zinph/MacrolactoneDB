import re
import sqlite3
# from flask import g
import pickle
import pandas as pd
# from flask_sqlalchemy import SQLAlchemy
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import IPythonConsole
from IPython.core.display import display, HTML
from flask import Flask, render_template, request, send_file, Response, url_for
from MacrolactoneDB_Miner import MacrolactoneDB_Miner

app = Flask(__name__)
# app.config['SQLALCHEMY_DATABASE_URI'] = "sqlite:///filtered_macrolactones.db"
# db = SQLAlchemy(app)


# def image_formatter(im):
#     return f'<img src="data:image/jpeg;base64,{image_base64(im)}">'

def path_to_image_html(path):
    return '<img src="'+ path + '" width="60" >'

def is_int(input):
  try:
    num = int(input)
  except ValueError:
    return False
  return True


def cleanup(input_string):
    '''
    When it's a range (-), it will be sorted into a tuple as min and max.
    If only one value is provided, both min and max are that value.
    '''
    stripped = re.sub(r'\s+', '', input_string)
    if stripped == None or stripped =='':
        low,high = 'dc','dc'
    elif '-' in stripped:
        splitted = stripped.split('-')
        low = int(min(splitted))
        high = int(max(splitted))
    else:
        if is_int(stripped):
            low,high = int(stripped), int(stripped)
        else:
            low,high = stripped, stripped
    return low,high


@app.route('/')
def user():
    img_file = url_for('static',filename='images/logo1.png')
    # img_file = "https://i.pinimg.com/originals/12/02/2e/12022e51379c7933a38add6c30bc7045.png"
    return render_template('user.html',img_file=img_file)

@app.route('/result',methods = ['POST', 'GET'])
def result():
   if request.method == 'POST':
      result = request.form
      return render_template("result.html",result = result)

@app.route('/library', methods=['POST'])
def getvalue():
    command = dict()
    ringsize = cleanup(request.form['ringsize'])
    command.update( [ ('RS_min', ringsize[0]) , ('RS_max', ringsize[1])])

    nRings = cleanup(request.form['nRings'])
    command.update( [ ('nRing_min', nRings[0]) , ('nRing_max', nRings[1])])

    nG12Rings = cleanup(request.form['nG12Rings'])
    command.update( [ ('nG12Ring_min', nG12Rings[0]) , ('nG12Ring_max', nG12Rings[1])])

    nFusedRings = cleanup(request.form['nFusedRings'])
    command.update( [ ('nFRing_min', nFusedRings[0]) , ('nFRing_max', nFusedRings[1])])

    nAromaticRings = cleanup(request.form['nAromaticRings'])
    command.update( [ ('naRing_min', nAromaticRings[0]) , ('naRing_max', nAromaticRings[1])])

    nSugars = cleanup(request.form['nSugars'])
    command.update( [ ('nSugars_min', nSugars[0]) , ('nSugars_max', nSugars[1])])

    nCoreEsters = cleanup(request.form['nCoreEsters'])
    command.update( [ ('core_ester_min', nCoreEsters[0]) , ('core_ester_max', nCoreEsters[1])])

    MW = cleanup(request.form['MW'])
    command.update( [ ('MW_min', MW[0]) , ('MW_max', MW[1])])

    SlogP = cleanup(request.form['SlogP'])
    command.update( [ ('SlogP_min', SlogP[0]) , ('SlogP_max', SlogP[1])])

    Lipinski = request.form['Lipinski']
    command['Lipinski'] = Lipinski

    activity_reported = request.form['activity_reported']
    command['activity_reported'] = activity_reported

    sample = MacrolactoneDB_Miner(command)
    library_df = sample.compile_filters()
    library_df.to_pickle('temp.pickle')
    # library_df.to_sql(name='temp', con=db.engine, index=False)
    num_cpds = library_df.shape[0]

    # smiles_writer(library_df['smiles'].tolist())
    # library_df = frame_manage(library_df)

    return render_template('view.html',tables=[library_df.to_html(index=False)], titles=library_df.columns.values, value=num_cpds)

@app.route('/return-smiles/')
def return_files_smiles():
    df = pd.read_pickle('temp.pickle')
    f = open('temp.smiles','w')
    f.write('smiles\n')
    f.write('\n'.join(df['smiles'].tolist()))
    f.close()
    try:
        result = send_file('temp.smiles', as_attachment = True)
        return result
    except Exception as e:
        return str(e)

@app.route('/return-CSV/')
def return_files_csv():
    df = pd.read_pickle('temp.pickle')
    df = df.drop(columns=['structures'])
    df.to_csv('temp.csv', index=False)
    try:
        result = send_file('temp.csv', as_attachment = True)
        return result
    except Exception as e:
        return str(e)

@app.route('/return-SDF/')
def return_files_sdf():
    df = pd.read_pickle('temp.pickle')
    PandasTools.WriteSDF(df, 'temp.sdf', molColName='structures', properties=list(df.columns), allNumeric=False)
    try:
        result = send_file('temp.sdf', as_attachment = True)
        return result
    except Exception as e:
        return str(e)


if __name__ == '__main__':
    app.run(debug=True)
   # app.run(port=80,debug = True)
