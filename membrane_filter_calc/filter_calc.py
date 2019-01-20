import pandas as pd
import numpy as np


def e_test(df):
    """
    This function determines if e condition has been met
    :param df:
    :return:
    """
    if df['Coliforms'].unique()[0] == 'mFC':
        typ_range = [20, 80]
    else:
        typ_range = [20, 60]

    e_value = ~df['num_Typical_Colonies'].between(*typ_range) & (df['num_Atypical_Colonies'] < 200) & (df['num_Typical_Colonies'] > 0)
    return e_value


def shade(df):
    return_ = """
    This function determines the row to be shaded
    :param df: 
    :return: 
    """
    df_temp = df.sort_values(by=['Invalid', 'Dilution Factor', 'Volume mL'], ascending=[True, True, False])
    shade_index = df_temp.index[0]

    return shade_index


def cfu_calc(mEndo_typ, mEndo_atyp, mFc_typ, mFc_atyp):
    df_meas_mEndo = pd.DataFrame.from_dict({'Dilution Factor': [1, 1, 1, 0.1, 0.01, 0.001], 'Volume mL': [0.5, 5, 50, 0.5, 5, 50]})
    df_meas_mEndo['Coliforms'] = 'mEndo'

    df_meas_mFC = pd.DataFrame.from_dict({'Dilution Factor': [1, 1, 1, 0.1, 0.01, 0.001], 'Volume mL': [0.5, 5, 50, 0.5, 5, 50]})
    df_meas_mFC['Coliforms'] = 'mFC'

    df_meas = pd.concat([df_meas_mEndo, df_meas_mFC], ignore_index=True)

    df_meas['Typical Colonies'] = ''
    df_meas['Atypical Colonies'] = ''

    df_meas.loc[df_meas['Coliforms'] == 'mEndo', 'Typical Colonies'] = np.array(mEndo_typ).astype('str')
    df_meas.loc[df_meas['Coliforms'] == 'mEndo', 'Atypical Colonies'] = np.array(mEndo_atyp).astype('str')

    df_meas.loc[df_meas['Coliforms'] == 'mFC', 'Typical Colonies'] = np.array(mFc_typ).astype('str')
    df_meas.loc[df_meas['Coliforms'] == 'mFC', 'Atypical Colonies'] = np.array(mFc_atyp).astype('str')

    # Replace '>200' by 201
    df_meas['num_Typical_Colonies'] = np.where(df_meas['Typical Colonies'] == '>200', 201, df_meas['Typical Colonies']).astype('int')
    df_meas['num_Atypical_Colonies'] = np.where(df_meas['Atypical Colonies'] == '>200', 201, df_meas['Atypical Colonies']).astype('int')

    df_meas['Total Colonies'] = df_meas['num_Typical_Colonies'] + df_meas['num_Atypical_Colonies']

    # Invalid Calculation
    df_meas['Invalid'] = df_meas['Total Colonies'] > 200

    # Invalid overwrites all
    # e calculation
    df_meas['e'] = df_meas.groupby('Coliforms').apply(e_test).values & (~df_meas['Invalid'])

    # < calculation
    df_meas['less_than'] = (df_meas['num_Typical_Colonies'] == 0) & (~df_meas['Invalid'])

    df_meas['CFU/mL'] = df_meas['num_Typical_Colonies'] * df_meas['Dilution Factor'] / df_meas['Volume mL'] * 100
    df_meas['results'] = df_meas['CFU/mL'].astype('str')

    df_meas['results'] = df_meas['results'].where(~df_meas['e'], df_meas['results'] + 'e')
    df_meas['results'] = df_meas['results'].where(~df_meas['less_than'], '<' + df_meas['results'])
    df_meas['results'] = df_meas['results'].where(~df_meas['Invalid'], 'Invalid')

    df_meas['shade'] = False
    shade_index = df_meas.groupby('Coliforms').apply(shade)
    df_meas.loc[shade_index.values, 'shade'] = True

    df_meas = df_meas.drop(columns=['CFU/mL']).rename(columns={'results': 'CFU/mL'})

    return df_meas


if __name__ == '__main__':
    mEndo_typ = [0, 0, 4, 21, 0, 0]
    mEndo_atyp = [1, 7, 17, 0, 0, 0]

    mFc_typ = [0, 0, 1, 0, 0, 0]
    mFc_atyp = [0, 0, 0, 0, 0, 0]

    results = cfu_calc(mEndo_typ, mEndo_atyp, mFc_typ, mFc_atyp)
