#!/usr/bin/env python3
"""View 1/4 — SPECIES-level x per-sample RELATIVE ABUNDANCE (phage->host genus).
Output: result_260605/fig4/test/view_species_relabund/
Run: conda run -n shotgun_virome python result_260605/fig4/code/render_view_species_relabund.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from _fig4_persample_render import render_view
NX = "/home/share/programs/nexvirome"
render_view("species", "relabund",
            f"{NX}/result_260605/fig4/test/view_species_relabund",
            "species_relabund")
