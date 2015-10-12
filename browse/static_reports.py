"""Module to load generated motif report pages."""

import os
import pickle

from collectf import settings

def get_static_reports(filename):
    """Gets the reports and ensemble report."""
    print 'getting reports'
    reports = pickle.load(open(os.path.join(
        settings.PICKLE_ROOT, 'reports', filename + '.pkl')))
    ensemble_report = pickle.load(open(os.path.join(
        settings.PICKLE_ROOT, 'ensemble_reports', filename + '.pkl')))
    return reports, ensemble_report
