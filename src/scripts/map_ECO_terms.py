# Map ECO terms to CollecTF experimental techniques

import pandas as pd

from core import models


def run():
    df = pd.read_csv('scripts/eco_terms.csv')
    for i, row in df.iterrows():
        technique_name = row['CollecTF experimental technique']
        ECO_term = row['ECO term ID to be mapped by CollecTF']
        # Find exp technique object in CollecTF database
        try:
            exp = models.ExperimentalTechnique.objects.get(name=technique_name)
            exp.EO_term = ECO_term
            exp.save()
        except models.ExperimentalTechnique.DoesNotExist:
            print technique_name, '- not found in CollecTF'

    print "Following techniques do not have a ECO term"
    print models.ExperimentalTechnique.objects.filter(EO_term__isnull=True)
    print models.ExperimentalTechnique.objects.filter(EO_term='')
