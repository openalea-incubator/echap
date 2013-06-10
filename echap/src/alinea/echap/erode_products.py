def erode_products(erosion_rates, compound_parameters = products_parameters):
    """ 
    Compute loss of activity of compounds due to evolution of resistance of strains
    :Parameters:
      - `erosion_rates` (dict) - A dict of ('compound': erosion_rate) where:
          - `erosion_rate` (float, [0,1]) - Rate of loss of efficacy of compounds due to evolution of resistance of strains
      - `coumpound_parameters` :see 'global_efficiency'
          
    :Returns:
      - `new_pars` : Same dict as `coumpound_parameters` with updated `Ap` and `Ae` linked to erosion of efficacy due to evolution of resistance of strains
      
    """
    
    ptable = dict([(p['compound'],p) for p in compound_parameters])
    new_pars = ptable.copy()
    for k in erosion_rates.keys() :
        if not k in ptable.keys():
            raise SimcyclePesticide('%s not found in compound_parameters'%k)
        
        new_pars[k]['Ap'] *= (1 - erosion_rates[k])
        new_pars[k]['Ae'] *= (1 - erosion_rates[k])

    return ptable.values(),