def grepunit(s):
  if '_emi' in s:
    return('emi')
  elif '_GMD' in s:
    return('GMD')
  elif '_ln(s)' in s:
    return('GSD')
  else: return s
