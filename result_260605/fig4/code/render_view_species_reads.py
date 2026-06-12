#!/usr/bin/env python3
"""View 3/4 — SPECIES-level x per-sample READ COUNT (phage->host genus).
Output: result_260605/fig4/test/view_species_reads/
Run: conda run -n shotgun_virome python result_260605/fig4/code/render_view_species_reads.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from _fig4_persample_render import render_view
NX = "/home/share/programs/nexvirome"
render_view("species", "reads",
            f"{NX}/result_260605/fig4/test/view_species_reads",
            "species_reads")
