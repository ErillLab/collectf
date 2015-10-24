"""Module to load generated motif report pages."""

import os
import pickle

from collectf import settings

def get_static_reports(filename):
    """Gets the reports and ensemble report."""
    report_data = pickle.load(open(os.path.join(
        settings.PICKLE_ROOT, 'reports', filename + '.pkl')))
    return report_data
