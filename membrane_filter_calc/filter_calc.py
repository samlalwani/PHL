import pandas as pd
import numpy as np
import ipywidgets as widgets
from ipywidgets import Layout
from IPython.display import display


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
    """
    This function determines the row to be shaded
    :param df: 
    :return: 
    """
    df_temp = df.sort_values(by=['Invalid', 'less_than', 'e', 'Dilution Factor', 'Volume mL'], ascending=[True, True, True, False, False])
    shade_index = df_temp.index[0]

    return shade_index


def cfu_calc(mEndo_typ, mEndo_atyp, mFc_typ, mFc_atyp):
    df_meas_mEndo = pd.DataFrame.from_dict(
        {'Dilution Factor': [1, 1, 1, 0.1, 0.01, 0.001], 'Volume mL': [0.5, 5, 50, 0.5, 5, 50], 'Typical Colonies': np.array(mEndo_typ).astype('str'), 'Atypical Colonies': np.array(mEndo_atyp).astype('str')})
    df_meas_mEndo['Coliforms'] = 'mEndo'

    df_meas_mFC = pd.DataFrame.from_dict(
        {'Dilution Factor': [1, 1, 1, 0.1, 0.01, 0.001], 'Volume mL': [0.5, 5, 50, 0.5, 5, 50], 'Typical Colonies': np.array(mFc_typ).astype('str'), 'Atypical Colonies': np.array(mFc_atyp).astype('str')})
    df_meas_mFC['Coliforms'] = 'mFC'

    df_meas = pd.concat([df_meas_mEndo, df_meas_mFC], ignore_index=True)

    # Replace '>200' by 201
    df_meas['num_Typical_Colonies'] = np.where(df_meas['Typical Colonies'] == '>200', 201, df_meas['Typical Colonies']).astype('int')
    df_meas['num_Atypical_Colonies'] = np.where(df_meas['Atypical Colonies'] == '>200', 201, df_meas['Atypical Colonies']).astype('int')

    df_meas['Total Colonies'] = df_meas['num_Typical_Colonies'] + df_meas['num_Atypical_Colonies']

    # Reporting Limit
    # TODO Check equation
    df_meas['Reporting Limit'] = 100 / df_meas['Volume mL'] / df_meas['Dilution Factor']

    # Invalid Calculation
    df_meas['Invalid'] = df_meas['Total Colonies'] > 200

    # Invalid overwrites all
    # e calculation
    df_meas['e'] = df_meas.groupby('Coliforms').apply(e_test).values & (~df_meas['Invalid'])

    # < calculation
    df_meas['less_than'] = (df_meas['num_Typical_Colonies'] == 0) & (~df_meas['Invalid'])

    # TODO Check calculation in documentation
    df_meas['CFU/mL'] = df_meas['num_Typical_Colonies'] * (100 / df_meas['Volume mL']) / df_meas['Dilution Factor']
    df_meas['results'] = df_meas['CFU/mL'].astype('str')

    df_meas['results'] = df_meas['results'].where(~df_meas['e'], df_meas['results'] + 'e')
    df_meas['results'] = df_meas['results'].where(~df_meas['less_than'], '<' + df_meas['results'])
    df_meas['results'] = df_meas['results'].where(~df_meas['Invalid'], 'Invalid')

    df_meas['shade'] = False
    shade_index = df_meas.groupby('Coliforms').apply(shade)
    df_meas.loc[shade_index.values, 'shade'] = True

    df_meas = df_meas.drop(columns=['CFU/mL']).rename(columns={'results': 'CFU/mL'})

    return df_meas


# The following functions are used to read data from widgets & to display the results

def highlight_row(row):
    '''
    Highlight the row where shade = True.
    '''
    if row['shade']:
        return ['background-color: yellow' for _ in range(row.shape[0])]
    else:
        return ['background-color: white' for _ in range(row.shape[0])]


def create_input_widget():
    """
    This function creates widget to read the inputs
    :return:
    """
    box_dilution_factor = widgets.VBox([widgets.Label(value="Dilution Factor", layout=Layout(width='15'))] + [widgets.Label(value=str(i)) for i in [1, 1, 1, 0.1, 0.01, 0.001]])
    box_volume_ml = widgets.VBox([widgets.Label(value="Volume mL")] + [widgets.Label(value=str(i)) for i in [0.5, 5, 50]] * 2)

    box_mEndo_typ = widgets.VBox([widgets.Label(value="mEndo Typical")] + [widgets.Text(value='0', layout=Layout(width='auto')) for i in range(6)])
    box_mEndo_atyp = widgets.VBox([widgets.Label(value="mEndo Atypical")] + [widgets.Text(value='0', layout=Layout(width='auto')) for i in range(6)])
    box_mFC_typ = widgets.VBox([widgets.Label(value="mFC Typical")] + [widgets.Text(value='0', layout=Layout(width='auto')) for i in range(6)])
    box_mFC_atyp = widgets.VBox([widgets.Label(value="mFC Atypical")] + [widgets.Text(value='0', layout=Layout(width='auto')) for i in range(6)])

    input_widget = widgets.HBox([box_dilution_factor, box_volume_ml, box_mEndo_typ, box_mEndo_atyp, box_mFC_typ, box_mFC_atyp], border=True)

    return input_widget


def display_results(input_widget):
    """

    :param input_widget:
    :return:
    """
    col = 2
    col_val = input_widget.children[col]
    mEndo_typ = [item.value for item in list(col_val.children)[1:]]

    col = 3
    col_val = input_widget.children[col]
    mEndo_atyp = [item.value for item in list(col_val.children)[1:]]

    col = 4
    col_val = input_widget.children[col]
    mFc_typ = [item.value for item in list(col_val.children)[1:]]

    col = 5
    col_val = input_widget.children[col]
    mFc_atyp = [item.value for item in list(col_val.children)[1:]]

    results = cfu_calc(mEndo_typ, mEndo_atyp, mFc_typ, mFc_atyp)

    print("mEndo Inputs")
    summary_mEndo = results.loc[results.Coliforms == 'mEndo', ['Dilution Factor', 'Volume mL', 'Typical Colonies', 'Atypical Colonies', 'shade']]
    summary_mEndo = summary_mEndo.style.apply(highlight_row, axis=1).hide_index()
    results_mEndo = results.loc[(results.Coliforms == 'mEndo') & results.shade, ['Dilution Factor', 'Volume mL', 'Typical Colonies', 'Atypical Colonies', 'Total Colonies', 'CFU/mL', 'Reporting Limit']].style.hide_index()

    summary_mFC = results.loc[results.Coliforms == 'mFC', ['Dilution Factor', 'Volume mL', 'Typical Colonies', 'Atypical Colonies', 'shade']]
    summary_mFC = summary_mFC.style.apply(highlight_row, axis=1).hide_index()
    results_mFC = results.loc[(results.Coliforms == 'mFC') & results.shade, ['Dilution Factor', 'Volume mL', 'Typical Colonies', 'Atypical Colonies', 'Total Colonies', 'CFU/mL', 'Reporting Limit']].style.hide_index()

    display(summary_mEndo)
    display(results_mEndo)

    print("mFC Inputs")
    display(summary_mFC)
    display(results_mFC)


if __name__ == '__main__':
    # Case 1
    mEndo_typ = [0, 0, 4, 21, 0, 0]
    mEndo_atyp = [1, 7, 17, 0, 0, 0]

    mFc_typ = [0, 0, 1, 0, 0, 0]
    mFc_atyp = [0, 0, 0, 0, 0, 0]

    # Case 2
    mEndo_typ = [20, '>200', '>200', 5, 2, 0]
    mEndo_atyp = [46, '>200', '>200', 6, 1, 0]

    mFc_typ = [0, 11, 109, 1, 0, 0]
    mFc_atyp = [0, 0, 0, 0, 0, 0]

    results = cfu_calc(mEndo_typ, mEndo_atyp, mFc_typ, mFc_atyp)

    a = 1
