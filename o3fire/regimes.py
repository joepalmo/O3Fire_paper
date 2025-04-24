import numpy as np
import pandas as pd
from .baseline import split_clean

def regime_definition(df, enhancements_df, method='absolute_tracer', clean_method='percentile', percentile=0.4):
    """
    Define regimes based on tracer populations
    """
    possible_methods = ['absolute_tracer', 'relative_tracer', 'absolute_tracer_CO', 'relative_tracer_CO']
    if method not in possible_methods:
        raise ValueError("method must be one of {} ; not '{}'".format(possible_methods, method))
        
    regime_df = df.copy(deep=True)
    regime_df['regime'] = np.nan
    
    #clean vs polluted
    #split up clean data
    clean_df = split_clean(df, method=clean_method, percentile=percentile)
    regime_df.loc[clean_df.index, 'regime'] = 'clean'
    
    if method == 'absolute_tracer':
        polluted_condition = (regime_df['regime']!='clean')
        extreme_fire_condition = ((regime_df['CH3CN'] >= regime_df.loc[polluted_condition,'CH3CN'].quantile(0.97))&(regime_df['C2Cl4'] <= regime_df.loc[polluted_condition,'C2Cl4'].quantile(0.85)))
        urban_condition = ((regime_df['C2Cl4'] > regime_df.loc[polluted_condition,'C2Cl4'].quantile(0.5)))
        secondary_urban_condition = (regime_df['CH2Cl2'] > regime_df.loc[polluted_condition, 'CH2Cl2'].quantile(0.5))
        extreme_urban_condition = ((regime_df['CH3CN'] <= regime_df.loc[polluted_condition,'CH3CN'].quantile(0.85))&(regime_df['C2Cl4'] >= regime_df.loc[polluted_condition,'C2Cl4'].quantile(0.97)))
        fire_condition = ((regime_df['CH3CN'] > regime_df.loc[polluted_condition,'CH3CN'].quantile(0.5)))
        secondary_fire_condition = (regime_df['HCN'] > regime_df.loc[polluted_condition, 'HCN'].quantile(0.5))
        
        regime_df.loc[(polluted_condition & extreme_fire_condition)|(polluted_condition & fire_condition & ~urban_condition), 'regime'] = 'fire'
        regime_df.loc[(polluted_condition & extreme_urban_condition)|(polluted_condition & urban_condition & ~fire_condition), 'regime'] = 'urban'
        regime_df.loc[polluted_condition & fire_condition & urban_condition, 'regime'] = 'heavy_mixed'
        regime_df.loc[polluted_condition & ~fire_condition & ~urban_condition, 'regime'] = 'light_mixed'
        
        regime_df.loc[(polluted_condition)&(((regime_df['CH3CN'].isna()))|((regime_df['C2Cl4'].isna()))), 'regime'] = 'none'
        
        secondary_condition = (regime_df['regime']=='none')
        regime_df.loc[(polluted_condition & secondary_condition & secondary_fire_condition & ~secondary_urban_condition), 'regime'] = 'fire'
        regime_df.loc[(polluted_condition & secondary_condition & ~secondary_fire_condition & secondary_urban_condition), 'regime'] = 'urban'
        regime_df.loc[polluted_condition & secondary_condition & secondary_fire_condition & secondary_urban_condition, 'regime'] = 'heavy_mixed'
        regime_df.loc[polluted_condition & secondary_condition & ~secondary_fire_condition & ~secondary_urban_condition, 'regime'] = 'light_mixed'
        
        regime_df.loc[(polluted_condition)&(((regime_df['CH3CN'].isna())&(regime_df['HCN'].isna()))|((regime_df['C2Cl4'].isna())&(regime_df['CH2Cl2'].isna()))), 'regime'] = 'none'
        
    elif method == 'relative_tracer':
        polluted_condition = (regime_df['regime']!='clean')
        extreme_fire_condition = ((enhancements_df['CH3CN'] >= enhancements_df.loc[polluted_condition,'CH3CN'].quantile(0.97))&(enhancements_df['C2Cl4'] <= enhancements_df.loc[polluted_condition,'C2Cl4'].quantile(0.85)))
        urban_condition = ((enhancements_df['C2Cl4'] > enhancements_df.loc[polluted_condition,'C2Cl4'].quantile(0.5)))
        secondary_urban_condition = (enhancements_df['CH2Cl2'] > enhancements_df.loc[polluted_condition, 'CH2Cl2'].quantile(0.5))
        extreme_urban_condition = ((enhancements_df['CH3CN'] <= enhancements_df.loc[polluted_condition,'CH3CN'].quantile(0.85))&(enhancements_df['C2Cl4'] >= enhancements_df.loc[polluted_condition,'C2Cl4'].quantile(0.97)))
        fire_condition = ((enhancements_df['CH3CN'] > enhancements_df.loc[polluted_condition,'CH3CN'].quantile(0.5)))
        secondary_fire_condition = (enhancements_df['HCN'] > enhancements_df.loc[polluted_condition, 'HCN'].quantile(0.5))
        
        regime_df.loc[(polluted_condition & extreme_fire_condition)|(polluted_condition & fire_condition & ~urban_condition), 'regime'] = 'fire'
        regime_df.loc[(polluted_condition & extreme_urban_condition)|(polluted_condition & urban_condition & ~fire_condition), 'regime'] = 'urban'
        regime_df.loc[polluted_condition & fire_condition & urban_condition, 'regime'] = 'heavy_mixed'
        regime_df.loc[polluted_condition & ~fire_condition & ~urban_condition, 'regime'] = 'light_mixed'
        
        regime_df.loc[(polluted_condition)&(((enhancements_df['CH3CN'].isna()))|((enhancements_df['C2Cl4'].isna()))), 'regime'] = 'none'
        
        secondary_condition = (regime_df['regime']=='none')
        regime_df.loc[(polluted_condition & secondary_condition & secondary_fire_condition & ~secondary_urban_condition), 'regime'] = 'fire'
        regime_df.loc[(polluted_condition & secondary_condition & ~secondary_fire_condition & secondary_urban_condition), 'regime'] = 'urban'
        regime_df.loc[polluted_condition & secondary_condition & secondary_fire_condition & secondary_urban_condition, 'regime'] = 'heavy_mixed'
        regime_df.loc[polluted_condition & secondary_condition & ~secondary_fire_condition & ~secondary_urban_condition, 'regime'] = 'light_mixed'
        
        regime_df.loc[(polluted_condition)&(((enhancements_df['CH3CN'].isna())&(enhancements_df['HCN'].isna()))|((enhancements_df['C2Cl4'].isna())&(enhancements_df['CH2Cl2'].isna()))), 'regime'] = 'none'
     
    elif method == 'absolute_tracer_CO':
        regime_df['CH3CN/CO'] = regime_df['CH3CN']/regime_df['CO']
        regime_df['C2Cl4/CO'] = regime_df['C2Cl4']/regime_df['CO']
        regime_df['HCN/CO'] = regime_df['HCN']/regime_df['CO']
        regime_df['CH2Cl2/CO'] = regime_df['CH2Cl2']/regime_df['CO']
        
        polluted_condition = (regime_df['regime']!='clean')
        extreme_fire_condition = ((regime_df['CH3CN/CO'] >= regime_df.loc[polluted_condition,'CH3CN/CO'].quantile(0.97))&(regime_df['C2Cl4/CO'] <= regime_df.loc[polluted_condition,'C2Cl4/CO'].quantile(0.85)))
        urban_condition = ((regime_df['C2Cl4/CO'] > regime_df.loc[polluted_condition,'C2Cl4/CO'].quantile(0.5)))
        secondary_urban_condition = (regime_df['CH2Cl2/CO'] > regime_df.loc[polluted_condition, 'CH2Cl2/CO'].quantile(0.5))
        extreme_urban_condition = ((regime_df['CH3CN/CO'] <= regime_df.loc[polluted_condition,'CH3CN/CO'].quantile(0.85))&(regime_df['C2Cl4'] >= regime_df.loc[polluted_condition,'C2Cl4/CO'].quantile(0.97)))
        fire_condition = ((regime_df['CH3CN/CO'] > regime_df.loc[polluted_condition,'CH3CN/CO'].quantile(0.5)))
        secondary_fire_condition = (regime_df['HCN/CO'] > regime_df.loc[polluted_condition, 'HCN/CO'].quantile(0.5))
        
        regime_df.loc[(polluted_condition & extreme_fire_condition)|(polluted_condition & fire_condition & ~urban_condition), 'regime'] = 'fire'
        regime_df.loc[(polluted_condition & extreme_urban_condition)|(polluted_condition & urban_condition & ~fire_condition), 'regime'] = 'urban'
        regime_df.loc[polluted_condition & fire_condition & urban_condition, 'regime'] = 'heavy_mixed'
        regime_df.loc[polluted_condition & ~fire_condition & ~urban_condition, 'regime'] = 'light_mixed'
        
        regime_df.loc[(polluted_condition)&(((regime_df['CH3CN/CO'].isna()))|((regime_df['C2Cl4/CO'].isna()))), 'regime'] = 'none'
        
        secondary_condition = (regime_df['regime']=='none')
        regime_df.loc[(polluted_condition & secondary_condition & secondary_fire_condition & ~secondary_urban_condition), 'regime'] = 'fire'
        regime_df.loc[(polluted_condition & secondary_condition & ~secondary_fire_condition & secondary_urban_condition), 'regime'] = 'urban'
        regime_df.loc[polluted_condition & secondary_condition & secondary_fire_condition & secondary_urban_condition, 'regime'] = 'heavy_mixed'
        regime_df.loc[polluted_condition & secondary_condition & ~secondary_fire_condition & ~secondary_urban_condition, 'regime'] = 'light_mixed'
        
        regime_df.loc[(polluted_condition)&(((regime_df['CH3CN/CO'].isna())&(regime_df['HCN/CO'].isna()))|((regime_df['C2Cl4/CO'].isna())&(regime_df['CH2Cl2/CO'].isna()))), 'regime'] = 'none'
        
    elif method == 'relative_tracer_CO':
        enhancements_df['CH3CN/CO'] = enhancements_df['CH3CN']/enhancements_df['CO']
        enhancements_df['C2Cl4/CO'] = enhancements_df['C2Cl4']/enhancements_df['CO']
        enhancements_df['HCN/CO'] = enhancements_df['HCN']/enhancements_df['CO']
        enhancements_df['CH2Cl2/CO'] = enhancements_df['CH2Cl2']/enhancements_df['CO']
        
        polluted_condition = (regime_df['regime']!='clean')
        extreme_fire_condition = ((enhancements_df['CH3CN/CO'] >= enhancements_df.loc[polluted_condition,'CH3CN/CO'].quantile(0.97))&(enhancements_df['C2Cl4/CO'] <= enhancements_df.loc[polluted_condition,'C2Cl4/CO'].quantile(0.85)))
        urban_condition = ((enhancements_df['C2Cl4/CO'] > enhancements_df.loc[polluted_condition,'C2Cl4/CO'].quantile(0.5)))
        secondary_urban_condition = (enhancements_df['CH2Cl2/CO'] > enhancements_df.loc[polluted_condition, 'CH2Cl2/CO'].quantile(0.5))
        extreme_urban_condition = ((enhancements_df['CH3CN/CO'] <= enhancements_df.loc[polluted_condition,'CH3CN/CO'].quantile(0.85))&(enhancements_df['C2Cl4'] >= enhancements_df.loc[polluted_condition,'C2Cl4/CO'].quantile(0.97)))
        fire_condition = ((enhancements_df['CH3CN/CO'] > enhancements_df.loc[polluted_condition,'CH3CN/CO'].quantile(0.5)))
        secondary_fire_condition = (enhancements_df['HCN/CO'] > enhancements_df.loc[polluted_condition, 'HCN/CO'].quantile(0.5))
        
        regime_df.loc[(polluted_condition & extreme_fire_condition)|(polluted_condition & fire_condition & ~urban_condition), 'regime'] = 'fire'
        regime_df.loc[(polluted_condition & extreme_urban_condition)|(polluted_condition & urban_condition & ~fire_condition), 'regime'] = 'urban'
        regime_df.loc[polluted_condition & fire_condition & urban_condition, 'regime'] = 'heavy_mixed'
        regime_df.loc[polluted_condition & ~fire_condition & ~urban_condition, 'regime'] = 'light_mixed'
        
        regime_df.loc[(polluted_condition)&(((enhancements_df['CH3CN/CO'].isna()))|((enhancements_df['C2Cl4/CO'].isna()))), 'regime'] = 'none'
        
        secondary_condition = (regime_df['regime']=='none')
        regime_df.loc[(polluted_condition & secondary_condition & secondary_fire_condition & ~secondary_urban_condition), 'regime'] = 'fire'
        regime_df.loc[(polluted_condition & secondary_condition & ~secondary_fire_condition & secondary_urban_condition), 'regime'] = 'urban'
        regime_df.loc[polluted_condition & secondary_condition & secondary_fire_condition & secondary_urban_condition, 'regime'] = 'heavy_mixed'
        regime_df.loc[polluted_condition & secondary_condition & ~secondary_fire_condition & ~secondary_urban_condition, 'regime'] = 'light_mixed'
        
        regime_df.loc[(polluted_condition)&(((enhancements_df['CH3CN/CO'].isna())&(enhancements_df['HCN/CO'].isna()))|((enhancements_df['C2Cl4/CO'].isna())&(enhancements_df['CH2Cl2/CO'].isna()))), 'regime'] = 'none'

    
    regime_df.loc[regime_df['regime'].isna(), 'regime'] = 'none'
        
    return regime_df['regime']