import pandas as pd

ui2ignore = [3473724, 3474352, 3481948, 3493671]

chal_dict = {
    "testing": {
        "bottlenack": 9616012,
        "secondary_contact": 9616013,
        "growth": 9616014,
        "split_mig": 9616015,
        "admix": 9616016,
        "single_sweep": 9616070,
        "single_sweep_bgs": 9616075,
        "multi_sweep": 9616073,
        "multi_sweep_bgs": 9616077,
    },
    "final": {
        "bottlenack": 9616018,
        "secondary_contact": 9616041,
        "growth": 9616042,
        "split_mig": 9616043,
        "admix": 9616044,
        "single_sweep": 9616072,
        "single_sweep_bgs": 9616076,
        "multi_sweep": 9616074,
        "multi_sweep_bgs": 9616078,
    }
}

chal_types = ["testing", "final"]
chals = ["bottlenack" ,"secondary_contact" ,"growth" ,"split_mig", "admix", "single_sweep", "single_sweep_bgs", "multi_sweep", "multi_sweep_bgs"]

chal_df = pd.read_csv("data/all-challenges-info.tsv", delimiter='\t', header=0)

# Remove uids of admin submissions
# ~ is for "not" to revearse true and false
chal_df = chal_df[~chal_df["createdBy"].isin(ui2ignore)]

for chal_type in chal_types:
    for chal in chals:
        chalid = chal_dict[chal_type][chal]
        # chal_type, chal, number of submissions, number of users that submitted
        print(chal_type,chal,len(chal_df[chal_df["evaluationid"]==chalid]),len(set(chal_df[chal_df["evaluationid"]==chalid]['createdBy'])))