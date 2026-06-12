#!/usr/bin/env python3
"""View 2/4 — GENUS-level x per-sample RELATIVE ABUNDANCE (phage->host genus).
Output: result_260605/fig4/test/view_genus_relabund/
Run: conda run -n shotgun_virome python result_260605/fig4/code/render_view_genus_relabund.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from _fig4_persample_render import render_view
NX = "/home/share/programs/nexvirome"
render_view("genus", "relabund",
            f"{NX}/result_260605/fig4/test/view_genus_relabund",
            "genus_relabund")
