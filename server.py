#!/usr/bin/env python3
"""
HLA-T2D Association Dashboard Server
======================================
Flask server for the 3D HLA-T2D visualization dashboard.

Usage:
    python server.py
    # Open http://localhost:5060 in your browser
"""

import json
import os
import tempfile
from pathlib import Path
from dataclasses import asdict

from flask import Flask, jsonify, request, send_file, send_from_directory
from flask_cors import CORS
from werkzeug.utils import secure_filename

from analyze_hla import (
    load_database,
    parse_alleles_string,
    parse_csv_input,
    generate_demo_typing,
    analyze_hla,
    DB_PATH,
)

app = Flask(__name__, static_folder="static")
CORS(app)

UPLOAD_FOLDER = tempfile.mkdtemp()
app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER
app.config["MAX_CONTENT_LENGTH"] = 50 * 1024 * 1024

db = load_database()


@app.route("/")
def index():
    return send_file("index.html")


@app.route("/static/<path:filename>")
def static_files(filename):
    return send_from_directory("static", filename)


@app.route("/api/database")
def get_database():
    return jsonify(db)


@app.route("/api/alleles")
def get_alleles():
    gene = request.args.get("gene")
    assoc = request.args.get("association")
    alleles = db["allele_associations"]
    if gene:
        alleles = [a for a in alleles if a["gene"] == gene]
    if assoc:
        alleles = [a for a in alleles if a["association"] == assoc]
    return jsonify({"alleles": alleles, "count": len(alleles)})


@app.route("/api/snps")
def get_snps():
    return jsonify({"snps": db["snp_associations"]})


@app.route("/api/haplotypes")
def get_haplotypes():
    return jsonify({"haplotypes": db["haplotype_risks"]})


@app.route("/api/mechanisms")
def get_mechanisms():
    return jsonify({"mechanisms": db["biological_mechanisms"]})


@app.route("/api/genes")
def get_genes():
    return jsonify({"genes": db["hla_genes"]})


@app.route("/api/analyze/demo", methods=["POST"])
def analyze_demo():
    typing = generate_demo_typing()
    report = analyze_hla(typing, db)
    return jsonify(asdict(report))


@app.route("/api/analyze/alleles", methods=["POST"])
def analyze_alleles():
    data = request.get_json()
    if not data or "alleles" not in data:
        return jsonify({"error": "Provide JSON with 'alleles' string"}), 400
    typing = parse_alleles_string(data["alleles"])
    if "snps" in data:
        typing.snps = data["snps"]
    report = analyze_hla(typing, db)
    return jsonify(asdict(report))


@app.route("/api/analyze/upload", methods=["POST"])
def analyze_upload():
    if "file" not in request.files:
        return jsonify({"error": "No file"}), 400
    file = request.files["file"]
    if not file.filename:
        return jsonify({"error": "Empty filename"}), 400
    filename = secure_filename(file.filename)
    filepath = os.path.join(app.config["UPLOAD_FOLDER"], filename)
    file.save(filepath)
    try:
        typing = parse_csv_input(filepath)
        report = analyze_hla(typing, db)
        return jsonify(asdict(report))
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    finally:
        if os.path.exists(filepath):
            os.remove(filepath)


if __name__ == "__main__":
    print("=" * 55)
    print("  HLA-T2D Association 3D Dashboard")
    print("  Open http://localhost:5060 in your browser")
    print("=" * 55)
    app.run(host="0.0.0.0", port=5060, debug=True)
