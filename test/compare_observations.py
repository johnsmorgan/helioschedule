import numpy as np
from astropy.table import Table
from yaml import safe_load

conf = safe_load(open("test.yaml"))

# %%
t1 = Table.read("observations_june17.csv")
t2 = Table.read("observations.csv")

# %%
print("lengths")
len(t1), len(t2)

# %%
print(t1.colnames)
print(t2.colnames)

# %%
print("==")
for col in t2.colnames:
        print(col, end=' ')
        print(np.all(np.nan_to_num(t1[col].data)==np.nan_to_num(t2[col].data)))
        #if not np.all(np.nan_to_num(t1[col].data)==np.nan_to_num(t2[col].data)):
        #        print('***************************')
print()
print("isclose")
for col in t2.colnames:
    if col=='local_noon_str':
        continue
    print(col, end=' ')
    #print(np.nan_to_num(t1[col].data).dtype, end=' ')
    #print(np.nan_to_num(t2[col].data).dtype, end=' ')
    print(np.all(np.isclose(np.nan_to_num(t1[col].data), np.nan_to_num(t2[col].data))))
print()
print("isclose beam, starttime")
for col in t2.colnames:
    if not (col.startswith('beam') or col.startswith('starttime')):
        continue
    print(col, end=' ')
    #print(np.nan_to_num(t1[col].data).dtype, end=' ')
    #print(np.nan_to_num(t2[col].data).dtype, end=' ')
    print(np.all(np.isclose(np.nan_to_num(t1[col].data), np.nan_to_num(t2[col].data))))
    print(t1[col].data)
    print(t2[col].data)

print()
print("others, excluding unflagged")
for col in t2.colnames:
    if (col.startswith('beam') or col.startswith('starttime') or col.startswith("unflagged")):
        continue
    if col=='local_noon_str':
        continue
    print(col, end=' ')
    #print(np.nan_to_num(t1[col].data).dtype, end=' ')
    #print(np.nan_to_num(t2[col].data).dtype, end=' ')
    a = np.all(np.isclose(np.nan_to_num(t1[col].data), np.nan_to_num(t2[col].data)))
    print(a)
    if not a:
        print(t1[col].data)
        print(t2[col].data)

print()
print("unflagged")
for col in t2.colnames:
    if not col.startswith("unflagged"):
        continue
    if col=='local_noon_str':
        continue
    print(col, end=' ')
    #print(np.nan_to_num(t1[col].data).dtype, end=' ')
    #print(np.nan_to_num(t2[col].data).dtype, end=' ')
    a = np.all(np.isclose(np.nan_to_num(t1[col].data), np.nan_to_num(t2[col].data)))
    print(a)
    if not a:
        print(t1[col].data)
        print(t2[col].data)
