/* Chemical LR Resource - static browse app.
   Same construction as the peptide site's app.js: no framework, no build step,
   data arrives as globals from data/*.js. Provenance is sharded and injected on
   demand when a row is opened. */
(function () {
  "use strict";
  const I = window.INTERACTIONS, C = window.CHEMICALS, E = window.ENZYMES,
        S = window.STATS, SRC = window.SRC_ORDER, SHARDS = S.shards;
  const PAGE = 100;
  const esc = s => String(s == null ? "" : s).replace(/[&<>"']/g,
    m => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[m]));

  /* Columns of an interaction row, kept in sync with np7_build_site.py:
     0 id  1 chemical_id  2 chemical_name  3 receptor  4 receptor_type
     5 n_sources  6 source index string  7 n_pmids  8 complex_id            */
  const SRC_LABEL = {
    gtopdb_interactions: "GtoPdb",
    cellchatdb_human: "CellChatDB v2",
    mrclinkdb_human_mr: "MRCLinkdb",
    cellphonedb_interactions: "CellPhoneDB v5",
    metaligand_pairs_human: "MetaLigand",
    metalinksdb: "MetalinksDB",
    metaligand_ligands_human: "MetaLigand", mrclinkdb_human_enzyme: "MRCLinkdb", mrclinkdb_human_transporter: "MRCLinkdb",
    human_gem: "Human-GEM v2.1.1", rhea_uniprot: "Rhea/UniProt"
  };
  /* Sources that carry a per-pair citation. Used only to colour the provenance
     strip - it is a provenance cue, never a score. */
  const CITING = new Set(["gtopdb_interactions", "cellchatdb_human", "mrclinkdb_human_mr"]);
  const ROLE_CLASS = { synthesis: "syn", degradation: "deg", efflux: "eff",
                       influx: "inf", precursor_synthesis: "pre" };
  const ROLE_HELP = {
    synthesis: "Genes whose products make the ligand itself.",
    precursor_synthesis: "Genes whose products make the ligand's precursor.",
    degradation: "Genes whose products consume or inactivate it.",
    efflux: "Transporters that export it, including synaptic-vesicle loading.",
    influx: "Transporters that import it, or its precursor (salvage)."
  };

  const srcList = bits => [...bits].map(i => SRC[+i]).filter(Boolean);
  const hasEnz = cid => !!(E[cid] && E[cid].length);

  function srcBar(bits) {
    const on = new Set(srcList(bits));
    return `<span class="srcbar" title="${esc(srcList(bits).map(s => SRC_LABEL[s] || s).join(", "))}">` +
      SRC.map(s => `<i class="${on.has(s) ? "on" + (CITING.has(s) ? " q" : "") : ""}"></i>`).join("") +
      `</span>`;
  }

  /* ---------------- complexes as single entries ----------------
     Decision 2026-10-19: a complex that the SOURCE databases mark (complex_id) is one
     entry per chemical, not one row per subunit. Only source-marked complexes are merged;
     our own family groupings are NOT (they were never reviewed). The merge happens here,
     at load time, so the exported CSVs and the database keep one row per subunit.
     A merged row keeps the shape of the others and adds  9 member ids, 10 member symbols. */
  (function collapseComplexes() {
    const groups = new Map();
    I.forEach(r => {
      if (!r[8]) return;
      const k = r[1] + "|" + r[8];
      if (!groups.has(k)) groups.set(k, []);
      groups.get(k).push(r);
    });
    if (!groups.size) return;
    const label = id => id.split(/[;_]/).filter(Boolean).join(" + ");
    const merge = (list, pick, keyf) => {
      const seen = new Set(), out = [];
      list.forEach(x => { const k = keyf(x); if (!seen.has(k)) { seen.add(k); out.push(x); } });
      return out;
    };
    const drop = new Set(), made = [];
    groups.forEach(rows => {
      rows.sort((a, b) => a[0] < b[0] ? -1 : 1);
      const first = rows[0], ids = rows.map(r => r[0]);
      const types = rows.map(r => r[4]).filter(Boolean);
      const m = [
        first[0], first[1], first[2], label(first[8]),
        types[0] || "",
        Math.max(...rows.map(r => r[5])),
        [...new Set(rows.flatMap(r => [...String(r[6])]))].sort().join(""),
        Math.max(...rows.map(r => r[7])),
        first[8], ids, rows.map(r => r[3])
      ];
      rows.forEach(r => drop.add(r));
      made.push(m);
      ["EVIDENCE", "CLASSEV"].forEach(name => {
        const store = window[name]; if (!store) return;
        const all = ids.flatMap(i => store[i] || []);
        if (all.length) store[first[0]] = merge(all, null, e => JSON.stringify([e.pmid, e.sent || e.det, e.cls]));
      });
    });
    const kept = I.filter(r => !drop.has(r));
    I.splice(0, I.length, ...kept, ...made);
  })();

  /* ---------- CSV download of the filtered table (owner request 2026-11-03) ---------- */
  function csvCell(v) {
    v = String(v == null ? "" : v);
    return /[",\n\r]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v;
  }
  function saveCSV(name, head, rows) {
    const text = "\ufeff" + [head, ...rows].map(r => r.map(csvCell).join(",")).join("\r\n");
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([text], { type: "text/csv;charset=utf-8" }));
    a.download = name; document.body.appendChild(a); a.click();
    setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
  }
  function downloadChemCSV() {
    const VS = new Set(window.VERIFIED || []);
    const head = ["pair_id", "chemical", "hmdb_id", "pubchem_cid", "formula", "chemical_class", "receptor", "receptor_type",
      "complex", "n_sources", "sources", "n_pmids", "verified", "literature_evidence", "assay_classes", "database_evidence", "n_genes"];
    saveCSV("lr_chemical_pairs.csv", head, view.map(r => {
      const ch = C[r[1]] || {}, ids = r[9] || [r[0]];
      const ev = (window.EVIDENCE || {})[r[0]] || [], cl = (window.CLASSEV || {})[r[0]] || [];
      return [r[0], r[2], ch.h || "", ch.pc || "", ch.f || "", ch.mc || "", r[3], r[4], r[8] || "", r[5],
        srcList(r[6]).map(x => SRC_LABEL[x] || x).join("; "), r[7], ids.some(i => VS.has(i)) ? "yes" : "no",
        [...new Set(ev.map(e => e.res).filter(Boolean))].join("; "), [...new Set(ev.map(e => e.as).filter(Boolean))].join("; "),
        [...new Set(cl.map(c => c.src).filter(Boolean))].join("; "), new Set((E[r[1]] || []).map(x => x[0])).size];
    }));
  }
  /* ---------------- filters + table ---------------- */
  let sortCol = 5, sortAsc = false, page = 0, view = I;

  function apply() {
    const q = document.getElementById("q").value.trim().toLowerCase();
    const mc = document.getElementById("mclass").value;
    const rt = document.getElementById("rtype").value;
    const ms = +document.getElementById("minsrc").value;
    const only = document.getElementById("only").value;
    const assay = document.getElementById("assay").value;
    view = I.filter(r => {
      if (q && !(r[2].toLowerCase().includes(q) || r[3].toLowerCase().includes(q) || (r[10] || []).some(g => g.toLowerCase().includes(q)))) return false;
      if (mc && (C[r[1]] || {}).mc !== mc) return false;
      if (rt && r[4] !== rt) return false;
      if (r[5] < ms) return false;
      if (only === "cited" && !r[7]) return false;
      if (only === "uncited" && r[7]) return false;
      if (only === "cplx" && !r[8]) return false;
      if (only === "enz" && !hasEnz(r[1])) return false;
      if (only === "noenz" && hasEnz(r[1])) return false;
      if (only === "prov" && !(C[r[1]] || {}).pv) return false;
      if (assay && !((window.EVIDENCE || {})[r[0]] || []).some(e => e.as === assay)) return false;
      if (only === "evid" && !(window.EVIDENCE || {})[r[0]]) return false;
      return true;
    });
    view.sort((a, b) => {
      let x = a[sortCol], y = b[sortCol];
      if (typeof x === "string") { x = x.toLowerCase(); y = y.toLowerCase(); }
      return (x < y ? -1 : x > y ? 1 : 0) * (sortAsc ? 1 : -1);
    });
    page = 0;
    render();
  }

  function render() {
    const n = view.length, pages = Math.max(1, Math.ceil(n / PAGE));
    if (page >= pages) page = pages - 1;
    const slice = view.slice(page * PAGE, page * PAGE + PAGE);
    document.getElementById("count").textContent =
      `${n.toLocaleString()} pair${n === 1 ? "" : "s"} of ${I.length.toLocaleString()} (a complex the sources mark counts as one)`;
    document.getElementById("rows").innerHTML = slice.map(r => {
      const ch = C[r[1]] || {}, ez = E[r[1]] || [];
      return `<tr data-id="${esc(r[0])}">
        <td class="pair">${esc(r[2])}${ch.pv ? ' <span class="chip" title="HMDB match may not be the exact stereoisomer or tautomer">provisional</span>' : ""}</td>
        <td class="pair">${esc(r[3])}${r[8] ? ` <span class="chip" title="subunit of ${esc(r[8])}">complex</span>` : ""}</td>
        <td>${r[4] ? `<span class="chip">${esc(r[4])}</span>` : ""}</td>
        <td>${srcBar(r[6])} <span style="color:var(--muted);font-size:12px">${r[5]}</span></td>
        <td>${r[7] ? esc(r[7]) : '<span class="unassessed" title="No upstream source cites a paper for this pair, and no quoted paper or database record carries a PMID">none</span>'}</td>
        <td>${evCell(r[0])}</td>
        <td>${assayCell(r[0])}</td>
        <td>${ch.mc ? `<span class="chip">${esc(ch.mc)}</span>` : ""}</td>
        <td>${ez.length ? esc(new Set(ez.map(x => x[0])).size) : '<span class="unassessed" title="No genes recorded for this chemical, so it cannot be inferred from expression data">none</span>'}</td>
      </tr>`;
    }).join("");
    document.getElementById("pageinfo").textContent = `page ${page + 1} / ${pages}`;
    document.getElementById("prev").disabled = page === 0;
    document.getElementById("next").disabled = page >= pages - 1;
  }

  /* ---------------- detail drawer ---------------- */
  const shardOf = id => { let s = 0; for (const c of id) s += c.charCodeAt(0); return s % SHARDS; };
  const loaded = {};
  function needShard(id, cb) {
    const ids = Array.isArray(id) ? id : [id];
    const todo = [...new Set(ids.map(shardOf))].filter(s => !loaded[s]);
    if (!todo.length) return cb();
    let left = todo.length;
    todo.forEach(s => {
      const el = document.createElement("script");
      el.src = `data/prov/p${s}.js`;
      el.onload = el.onerror = () => { loaded[s] = 1; if (--left === 0) cb(); };
      document.head.appendChild(el);
    });
  }


  /* Per-pair evidence summary in the pairs table, the same shape the peptide table
     uses: how many primary papers back this pair and what kinds of experiment they
     are. The quoted sentences stay in the drawer - a table cell cannot hold them. */
  const RES_SHORT = { SUPPORTS: "supports", REFUTES: "refutes", FAMILY_LEVEL: "family",
                      TRANSPORTER: "transport", ENZYME: "enzyme" };
  const ASSAY_SHORT = { structural: "struct", direct_binding: "binding",
                        functional_in_vitro: "in vitro", functional_in_vivo: "in vivo",
                        computational: "in silico" };
  /* Non-literature evidence classes (CITATION_METHODS.md Methods 3, 4). A
     nomenclature_evident pair ranks ABOVE literature evidence (owner decision
     2026-10-15); database_annotated ranks below it. Since 2026-10-31 (owner) a database_annotated
     record also makes the pair verified. */
  const CLS_LABEL = { nomenclature_evident: "name-evident", database_annotated: "database" };
  const CLS_TITLE = {
    nomenclature_evident: "The receptor's official HGNC name contains this molecule's own name (exact name only, no class matching). High confidence without a citation.",
    database_annotated: "A curated database (ChEMBL, GtoPdb, BindingDB, PDB, UniProt/Rhea or TCDB) records this molecule acting on, binding or being transported by this protein. Counted as direct evidence; no sentence has been quoted from a paper."
  };
  function clsBadges(id) {
    const cl = (window.CLASSEV || {})[id] || [];
    return [...new Set(cl.map(c => c.cls))].sort((a, b) => a === "nomenclature_evident" ? -1 : 1)
      .map(k => `<span class="badge b-verified" title="${esc(CLS_TITLE[k] || "")}">${esc(CLS_LABEL[k] || k)}</span>`).join("");
  }
  function classPanel(id) {
    const cl = (window.CLASSEV || {})[id] || [];
    if (!cl.length) return "";
    const order = c => c.cls === "nomenclature_evident" ? 0 : 1;
    return `<h3 style="margin:22px 0 8px;font-size:15px">Evidence class</h3>` +
      cl.slice().sort((a, b) => order(a) - order(b)).map(c => `<div class="ev confirmed">
        <div class="top">
          <span class="badge b-verified" title="${esc(CLS_TITLE[c.cls] || "")}">${esc(CLS_LABEL[c.cls] || c.cls)}</span>
          <span class="chip">${esc(c.conf)} confidence</span>
          <span class="chip">${esc(c.src)}</span>
          ${c.pmid ? `<span style="margin-left:auto">${pmidLinks(c.pmid)}</span>` : ""}
        </div>
        <div class="title" style="color:var(--muted);margin-top:6px">${esc(c.det)}</div>
      </div>`).join("");
  }
  function evCell(id) {
    const ev = (window.EVIDENCE || {})[id] || [];
    const cb = clsBadges(id);
    if (!ev.length) return cb || '<span class="unassessed" title="No primary paper has been read for this pair">none</span>';
    const res = [...new Set(ev.map(e => e.res).filter(Boolean))];
    const ref = res.includes("REFUTES");
    return cb + `<span class="badge b-verified" title="${ev.length} primary paper${ev.length === 1 ? "" : "s"} read; quoted sentence verified against the paper's own text">${ev.length} paper${ev.length === 1 ? "" : "s"}</span>` +
      res.map(k => `<span class="chip"${ref && k === "REFUTES" ? ' style="color:#8a2b2b;border-color:#e0a0a0"' : ""}>${esc(RES_SHORT[k] || k.toLowerCase())}</span>`).join("");
  }
  function assayCell(id) {
    const ev = (window.EVIDENCE || {})[id] || [];
    if (!ev.length) return "";
    const cl = [...new Set(ev.map(e => e.as).filter(Boolean))];
    return cl.map(c => `<span class="chip" title="${esc(c.replace(/_/g, " "))}">${esc(ASSAY_SHORT[c] || c)}</span>`).join("");
  }

  function pmidLinks(ref) {
    const clean = /^\s*\d+(\s*[;|,]\s*\d+)*\s*$/.test(ref || "");   // a plain id list: short old PMIDs (e.g. 8699) are real
    const ids = [...new Set((ref || "").match(clean ? /\d+/g : /(?<!\d)\d{7,8}(?!\d)/g) || [])];
    if (!ids.length) return "";
    return ids.map(p => `<a class="mono" href="https://pubmed.ncbi.nlm.nih.gov/${p}/"
      target="_blank" rel="noopener">PMID ${p} &#8599;</a>`).join(" ");
  }

  const GT_LABEL = { enzyme: "enzyme", transporter: "transporter", channel: "channel", other: "other" };
  const GT_HELP = {
    enzyme: "Catalyses a reaction (HGNC lists an EC number, or the protein name ends in -ase).",
    transporter: "Moves molecules across a membrane (solute carriers, ABC transporters, pumps).",
    channel: "Ion channel subunit.",
    other: "None of the above in HGNC's records (for example binding proteins and regulators)."
  };
  /* Root-source check (np22): does Human-GEM / Rhea / TCDB / HMDB back this gene in this role? */
  const RSRC = { G: "Human-GEM", R: "Rhea/UniProt", T: "TCDB", H: "HMDB", X: "Reactome", O: "GO (experimental)" };
  const OPP = { efflux: "influx", influx: "efflux" };
  function rootMark(cid, g, role) {
    const rr = ((window.GENEROOT || {})[cid] || {})[g] || [];
    const gr = rr.filter(x => /[GR]/.test(x[1]));          /* Human-GEM / Rhea: may set direction */
    const go = rr.filter(x => /O/.test(x[1]));             /* GO: supports, never sets direction */
    const names = x => [...x[1]].map(k => RSRC[k]).join(", ");
    const mk = (sym, col, tip) => `<span class="rmk" style="color:${col};font-weight:700;margin-left:3px" title="${esc(tip)}">${sym}</span>`;
    if (role === "precursor_synthesis") return mk("–", "var(--muted)", "Precursor genes are not checked against root sources on the site");
    let hit = gr.find(x => x[0] === role && x[2]);
    if (hit) return mk("✓", "#2e7d32", "Root agrees: " + names(hit) + " has a reaction with this gene in this role");
    /* a directed reaction in the opposite direction outranks an undirected match */
    if (OPP[role] && gr.find(x => x[0] === OPP[role] && x[2])) return mk("⇄", "#c62828", "Root says the OPPOSITE direction: " + OPP[role]);
    hit = gr.find(x => x[0] === role);
    if (hit) return mk("≈", "#8d6e00", "Root agrees only through a reversible or undirected reaction (" + names(hit) + ")");
    if (OPP[role] && rr.find(x => x[0] === "transport")) return mk("≈", "#8d6e00", "TCDB lists this gene as a transporter of the molecule (no direction)");
    if (go.find(x => x[0] === role || (OPP[role] && x[0] === "transport"))) return mk("✓", "#2e7d32", "GO (experimental evidence) links this gene to the molecule in this role; GO does not set a direction");
    if (gr.length) return mk("≠", "#c62828", "Root lists this gene in another role: " + [...new Set(gr.map(x => x[0]))].join(", "));
    if (go.length) return mk("≠", "#c62828", "GO (experimental) links this gene to the molecule in another role: " + [...new Set(go.map(x => x[0]))].join(", "));
    if (rr.find(x => /X/.test(x[1]))) return mk("◇", "#1565c0", "Reactome: gene and molecule take part in the same curated reaction (role not stated)");
    if (rr.find(x => /H/.test(x[1]))) return mk("○", "var(--muted)", "UNSUPPORTED: only HMDB links this gene and molecule (no role given)");
    return mk("·", "var(--muted)", "UNSUPPORTED: no source (Human-GEM, Rhea, Reactome, GO, TCDB, HMDB) links this gene to the molecule");
  }
  function rootOnly(cid, byRole) {
    const R = (window.GENEROOT || {})[cid] || {}, GI = window.GENEINFO || {};
    const add = {};
    Object.entries(R).forEach(([g, rr]) => rr.forEach(([role, m, dir]) => {
      if (!["synthesis", "degradation", "efflux", "influx"].includes(role) || !/[GR]/.test(m) || !dir) return;
      if (byRole[role] && byRole[role][g]) return;
      (add[role] = add[role] || {})[g] = m;
    }));
    const n = Object.values(add).reduce((t, o) => t + Object.keys(o).length, 0);
    if (!n) return "";
    return `<details style="margin-top:8px"><summary style="font-size:12.5px">Root sources name ${n} more gene${n === 1 ? "" : "s"} that no aggregated resource lists</summary>` +
      ["synthesis", "efflux", "influx", "degradation"].filter(r => add[r]).map(role =>
        `<div style="margin:6px 0"><div style="font-size:12px;color:var(--muted)">${esc(role)} · ${Object.keys(add[role]).length}</div>` +
        Object.entries(add[role]).sort().map(([g, m]) =>
          `<span class="pill" style="opacity:.75" title="${esc((GI[g] || ["", ""])[1] + " — from " + [...m].map(k => RSRC[k]).join(", "))}">${esc(g)}</span>`).join("") + `</div>`).join("") +
      `</details>`;
  }

  function enzPanel(cid) {
    const rows = E[cid] || [], GI = window.GENEINFO || {};
    if (!rows.length && (window.HUBS || []).includes(cid)) return `<p style="font-size:12.5px;padding:6px 8px;border-left:3px solid #8d6e00;background:var(--chip)">
      <b>Cofactor hub: no genes listed.</b> This molecule is made or used by hundreds of enzymes across metabolism, so no gene
      set can stand in for where it is released as a signal. Its availability cannot be estimated from expression data.</p>`;
    if (!rows.length) return `<p style="color:var(--muted)">No genes recorded. This chemical
      <b>cannot be inferred from transcriptomic data</b> — there is no gene whose expression
      stands in for it.</p>` + rootOnly(cid, {});
    const byRole = {};
    rows.forEach(([g, role, src, dk, also, origin]) => {
      (byRole[role] = byRole[role] || {});
      (byRole[role][g] = byRole[role][g] || { srcs: new Set(), dk: dk, also: new Set(), fb: origin === "lr_fallback" });
      byRole[role][g].srcs.add(src);
      (also || "").split("|").filter(Boolean).forEach(x => byRole[role][g].also.add(SRC_LABEL[x] || x));
    });
    const anyFb = rows.some(r => r[5] === "lr_fallback");
    const order = ["precursor_synthesis", "synthesis", "efflux", "influx", "degradation"];
    const isHub = (window.HUBS || []).includes(cid);
    const nGenes = new Set(rows.map(r => r[0])).size;
    return (isHub ? `<p style="font-size:12.5px;margin:4px 0 6px;padding:6px 8px;border-left:3px solid #8d6e00;background:var(--chip)">
      <b>Cofactor hub.</b> This molecule is used or made by hundreds of enzymes across metabolism, so these genes say nothing about
      where it is released as a signal. Do not use them as an expression proxy for its availability.</p>
      <details><summary style="font-size:12.5px">Show all ${nGenes} genes</summary>` : "") +
      `<p style="color:var(--muted);font-size:12px;margin:0 0 4px">Within each role, genes are ordered by how many
      raw sources give them (the number on each gene: Human-GEM and/or Rhea); the tag says what kind of protein it is.
      Genes come from the raw sources themselves, the Human-GEM metabolic model and Rhea/UniProt reactions, using only
      reactions with a stated direction; hover a gene to see which ligand&ndash;receptor resources also list it.</p>` +
      (anyFb ? `<p style="font-size:12.5px;margin:4px 0 6px;padding:6px 8px;border-left:3px solid var(--muted);background:var(--chip)">
      No raw source gives genes for this chemical, so the genes listed by ligand&ndash;receptor resources are shown instead.
      They are <b>not confirmed</b> by Human-GEM or Rhea (faded).</p>` : "") +

      order.filter(r => byRole[r]).map(role => {
      const genes = Object.entries(byRole[role]).sort((a, b) => b[1].srcs.size - a[1].srcs.size || (a[0] < b[0] ? -1 : 1));
      return `<div style="margin:8px 0">
        <div style="font-size:12px;color:var(--muted)" title="${esc(ROLE_HELP[role] || "")}">
          ${esc(role.replace("_", " "))} · ${genes.length} gene${genes.length === 1 ? "" : "s"}</div>
        <div style="margin-top:4px">${genes.map(([g, v]) => {
          const gi = GI[g] || ["other", ""], t = gi[0];
          const rm = v.fb ? rootMark(cid, g, role) : "", weak = v.fb;
          return `<span class="pill ${ROLE_CLASS[role] || ""}"${weak ? ' style="opacity:.5"' : ""} title="${esc((gi[1] ? gi[1] + " — " : "") + (v.fb ? "listed by " : "from ") + [...v.srcs].map(s => SRC_LABEL[s] || s).join(", ") + (v.also.size && !v.fb ? "; also listed by " + [...v.also].join(", ") : ""))}${v.dk ? "" : " — direction not stated by the source"}">${esc(g)}${v.dk ? "" : "?"}<sup class="gcount">${v.srcs.size}</sup><span class="gtype gt-${t}" title="${esc(GT_HELP[t])}">${esc(GT_LABEL[t])}</span>${rm}</span>`;
        }).join("")}</div></div>`;
    }).join("") + rootOnly(cid, byRole) + (isHub ? "</details>" : "");
  }

  /* ---------------- confirmed evidence panel ----------------
     chemical_evidence.csv rows for this pair. These are NOT the upstream citations
     below: each one is a primary paper whose quoted sentence was checked verbatim
     against the paper's own abstract or publisher full text at load time. They were
     confirmed without a curator reading the paper, so the label says so. */
  const RES_TITLE = {
    SUPPORTS: "the paper tests this molecule at this receptor and finds an interaction",
    REFUTES: "the paper tests this pair and finds NO interaction - the pair is contradicted",
    FAMILY_LEVEL: "the paper implicates the receptor family, not the named subunit",
    TRANSPORTER: "the protein moves the molecule; that is not receptor binding",
    ENZYME: "the protein consumes or produces the molecule; that is not receptor binding"
  };
  function evidencePanel(id, row) {
    const ev = (window.EVIDENCE || {})[id] || [];
    if (!ev.length) return "";
    return `<h3 style="margin:22px 0 8px;font-size:15px">Confirmed evidence &mdash; ${ev.length} paper${ev.length === 1 ? "" : "s"}</h3>` +
      ev.map(e => `<div class="ev confirmed">
        <div class="top">
          <span class="badge b-verified" title="${esc(RES_TITLE[e.res] || "")}">${esc((e.res || "unlabelled").replace(/_/g, " "))}</span>
          <span class="chip">${esc((e.as || "").replace(/_/g, " "))}</span>
          ${e.dir && e.dir !== "not determined" ? `<span class="chip">${esc(e.dir)}</span>` : ""}
          ${e.dep === "full_text" ? '<span class="chip" title="Quote taken from publisher full text, not the abstract">full text</span>' : ""}
          <span style="margin-left:auto">${pmidLinks(e.pmid)}</span>
        </div>
        <blockquote style="margin:8px 0 6px;padding-left:11px;border-left:3px solid var(--ok,#1a7f37);font-size:13.5px;line-height:1.6">
          &ldquo;${esc(e.sent)}&rdquo;</blockquote>
        <div class="title" style="color:var(--muted)">
          ${esc(e.mol || (row && row[2]) || "")} at ${esc(e.rec || (row && row[3]) || "")}${e.jr ? " &middot; " + esc(e.jr) : ""}${e.yr ? " " + esc(e.yr) : ""}
        </div>
        ${e.aff ? `<div class="title" style="color:var(--muted)">${esc(e.aff)}</div>` : ""}
      </div>`).join("") +
      `<p style="color:var(--muted);font-size:12px;margin:8px 0 0">Quoted sentence verified verbatim
       against the paper's own text by the load gates; confirmed without a curator read.</p>`;
  }

  /* HMDB context for the molecule (not for the pair): where it is found, in what
     sample, in which cell compartment, and which diseases HMDB links it to. */
  const ctxChips = a => a.map(v => `<span class="chip">${esc(v)}</span>`).join(" ");
  function ctxList(label, a, tip) {
    if (!a || !a.length) return "";
    const shown = a.length > 12
      ? ctxChips(a.slice(0, 12)) + ` <details style="display:inline"><summary>+${a.length - 12} more</summary>${ctxChips(a.slice(12))}</details>`
      : ctxChips(a);
    return `<dt title="${esc(tip)}">${label}</dt><dd style="line-height:1.9">${shown}</dd>`;
  }
  const uniq = a => a && [...new Set(a)];
  const DRUGPW = /\bAction\b/i;   /* SMPDB names drug-action pathways "<Drug> Action Pathway" */
  function pathRows(x) {
    const p = x.p;
    if (!p || !p.length) return "";
    const own = p.filter(n => !DRUGPW.test(n)), drug = p.filter(n => DRUGPW.test(n));
    const more = (a, label) => a.length
      ? ` <details style="display:inline"><summary>${label}</summary>${ctxChips(a)}</details>` : "";
    const head = own.slice(0, 10), rest = own.slice(10);
    return `<dt title="Pathways in HMDB / SMPDB that list this metabolite as a substrate, product, cofactor or other participant. A pathway listing does not mean the molecule signals through it.">Pathways</dt>
      <dd style="line-height:1.9">${head.length ? ctxChips(head) : '<span style="color:var(--muted);font-size:12.5px">none besides drug-action pathways</span>'}
        ${more(rest, "+" + rest.length + " more")}
        ${more(drug, drug.length + " drug-action pathway" + (drug.length === 1 ? "" : "s"))}
        <div style="color:var(--muted);font-size:12px">HMDB / SMPDB. Lipid-species variants of one pathway are merged.</div></dd>`;
  }
  function roleRows(x) {
    const l = x.l;
    if (!l) return "";
    return Object.keys(l).sort().map(g => ctxList(g, uniq(l[g]),
      "HMDB ontology: " + g)).join("");
  }
  function ctxRows(ch) {
    const x = ch.x;
    if (!x) return `<dt>Found in</dt><dd style="color:var(--muted);font-size:12.5px">No tissue, location or disease annotation in HMDB for this chemical.</dd>`;
    return ctxList("Tissues", x.t, "Tissues in which HMDB reports this metabolite") +
      ctxList("Biospecimens", x.b, "Body fluids or samples in which HMDB reports it") +
      ctxList("Cell location", x.c, "Cellular compartments HMDB lists for it") +
      ctxList("Diseases", x.d, "Diseases HMDB links to altered levels of this metabolite") +
      pathRows(x) +
      ctxList("Origin", uniq(x.s), "Where HMDB says the molecule comes from: made by the body (endogenous), food, microbes, plants, synthetic and so on") +
      ctxList("Route of exposure", uniq(x.r), "How HMDB says the body is exposed to it") +
      roleRows(x) +
      `<dt></dt><dd style="color:var(--muted);font-size:12px">From HMDB 5.0. These describe the molecule, not this ligand&ndash;receptor pair, and are not evidence for it.</dd>`;
  }

  function openDetail(id) {
    const row0 = I.find(x => x[0] === id);
    const members = (row0 && row0[9]) || [id];
    needShard(members, () => {
      const r = row0;
      if (!r) return;
      const ch = C[r[1]] || {};
      const prov = members.flatMap(m => (window.__PROV[shardOf(m)] || {})[m] || []);
      /* One line: which databases report this pair, and the PMIDs they cite between
         them. The per-source rows with the upstream spelling stay underneath in a
         disclosure, because the Methods tab promises that traceability. */
      const provNames = [...new Set(prov.map(p => SRC_LABEL[p[0]] || p[0]))];
      const provPmids = [...new Set(prov.flatMap(p =>
        (p[3] || "").match(/(?<!\d)\d{7,8}(?!\d)/g) || []))];
      const nCite = prov.filter(p => p[3]).length;
      const dbNames = [...new Set(((window.CLASSEV || {})[id] || []).map(c => c.src).filter(Boolean))];
      const allPm = [...new Set([...provPmids,
        ...((window.EVIDENCE || {})[id] || []).flatMap(e => String(e.pmid || "").match(/\d{1,9}/g) || []),
        ...((window.CLASSEV || {})[id] || []).flatMap(c => String(c.pmid || "").match(/\d{1,9}/g) || [])])];
      const provHtml = prov.length ? `
        <p style="margin:0 0 8px;font-size:13px">
          <b>${provNames.length} database${provNames.length === 1 ? "" : "s"}:</b>
          ${esc(provNames.join(" · "))}
          <span style="color:var(--muted)">&mdash; ${nCite
            ? nCite + " of " + prov.length + " row" + (prov.length === 1 ? "" : "s") + " cite a paper"
            : "none of these source rows cites a paper"}</span>
          ${dbNames.length ? `<br><b>Database evidence:</b> ${esc(dbNames.join(" · "))}` : ""}
          ${allPm.length ? "<br>" + allPm.map(x => pmidLinks(x)).join(" ") : ""}
        </p>
        <details style="font-size:12.5px">
          <summary style="cursor:pointer;color:var(--muted)">Upstream rows, as each database wrote them</summary>
          ${prov.map(p => `<div class="ev ${CITING.has(p[0]) ? "confirmed" : "pending"}" style="margin-top:6px">
            <div class="top">
              <span class="badge ${CITING.has(p[0]) ? "b-verified" : "b-unverified"}">${esc(SRC_LABEL[p[0]] || p[0])}</span>
              ${p[3] ? "" : '<span class="chip" title="This source supplies no citation for this pair">no citation</span>'}
              <span style="margin-left:auto">${pmidLinks(p[3])}</span>
            </div>
            <div class="title"><code>${esc(p[1])}</code> &rarr; <code>${esc(p[2])}</code></div>
          </div>`).join("")}
        </details>`
        : `<p style="color:var(--muted)">No source rows found.</p>`;

      const ext = [];
      if (ch.h) ext.push(`<a href="https://hmdb.ca/metabolites/${esc(ch.h)}" target="_blank" rel="noopener">${esc(ch.h)} &#8599;</a>`);
      if (ch.pc) ext.push(`<a href="https://pubchem.ncbi.nlm.nih.gov/compound/${esc(ch.pc)}" target="_blank" rel="noopener">PubChem ${esc(ch.pc)} &#8599;</a>`);

      document.getElementById("detail").innerHTML = `
      <div class="scrim"><div class="drawer">
        <button class="close">Close</button>
        <h2>${esc(r[2])} <span class="arrow">&rarr;</span> ${esc(r[3])}</h2>
        <div style="color:var(--muted);font-size:13px">
          ${esc(r[4] || "receptor type unrecorded")}${r[8] ? (r[9]
            ? ` · complex <code>${esc(r[8])}</code> shown as one entry; subunit rows merged: ${esc(r[10].join(", "))}`
            : ` · subunit of <code>${esc(r[8])}</code>`) : ""}
        </div>
        ${ch.pv ? `<div class="warn" style="margin-top:12px">The HMDB entry matched for this
          chemical may not be the exact stereoisomer or tautomer the source intended.</div>` : ""}

        <dl class="kv">
          <dt>Identifiers</dt><dd>${ext.join(" · ") || "<i>none — name-keyed only</i>"}
            ${ch.h ? "" : '<div style="color:var(--muted);font-size:12px">This chemical never resolved to an HMDB id; GtoPdb ships none.</div>'}</dd>
          <dt>ChEBI</dt><dd>&mdash;</dd>
          <dt>Formula</dt><dd>${esc(ch.f) || "&mdash;"}</dd>
          <dt>Class</dt><dd>${ch.mc ? `<span class="chip">${esc(ch.mc)}</span>` : "&mdash;"}
            ${ch.sc ? `<div style="color:var(--muted);font-size:12px">HMDB: ${esc(ch.sc)}${ch.cl ? " / " + esc(ch.cl) : ""}</div>` : ""}
            ${ch.cc ? `<div style="color:var(--muted);font-size:12px" title="MetaLigand's own grouping (KEGG-style). It lumps amino acids, amines and neurotransmitters under 'Peptides'; it is not a ChEBI class and not the chemical class used elsewhere on this site.">MetaLigand group: ${esc(ch.cc)}</div>` : ""}</dd>
          ${ch.sm ? `<dt title="SMILES is a one-line text code for a molecule's structure: letters are atoms, brackets and digits show branches and rings. Software can redraw the structure from it.">SMILES</dt><dd class="smiles">${esc(ch.sm)}</dd>` : ""}
          ${ctxRows(ch)}
          ${ch.syn && ch.syn.length ? `<dt>Also known as</dt><dd style="font-size:12.5px;color:var(--muted)">${esc(ch.syn.slice(0, 8).join(", "))}</dd>` : ""}
        </dl>

        <dl class="kv">${window.geneSection("Receptor gene links", [r[3]], "any")}</dl>

        <h3 style="margin:0 0 6px;font-size:15px">Related genes for the ligand (${esc(r[2])})</h3>
        ${enzPanel(r[1])}

        ${classPanel(id)}
        ${evidencePanel(id, r)}

        <h3 style="margin:22px 0 8px;font-size:15px">Reported by</h3>
        ${provHtml}
      </div></div>`;
      const close = () => { document.getElementById("detail").innerHTML = ""; };
      document.querySelector("#detail .close").onclick = close;
      document.querySelector("#detail .scrim").onclick = e => {
        if (e.target.classList.contains("scrim")) close();
      };
    });
  }

  /* ---------------- chemicals tab ---------------- */
  let cpage = 0, cview = [];
  const recCount = {};
  I.forEach(r => { recCount[r[1]] = (recCount[r[1]] || 0) + 1; });

  function capply() {
    const q = document.getElementById("cq").value.trim().toLowerCase();
    const only = document.getElementById("conly").value;
    cview = Object.entries(C).filter(([id, c]) => {
      if (q && !((c.n || "").toLowerCase().includes(q) ||
                 (c.h || "").toLowerCase().includes(q) ||
                 (c.f || "").toLowerCase().includes(q))) return false;
      if (only === "enz" && !hasEnz(id)) return false;
      if (only === "noenz" && hasEnz(id)) return false;
      if (only === "nt" && !c.nt) return false;
      if (only === "prov" && !c.pv) return false;
      const x = c.x || {};
      const ft = document.getElementById("ctissue").value,
            fb = document.getElementById("cbio").value,
            fd = document.getElementById("cdis").value.trim().toLowerCase();
      if (ft && !(x.t || []).includes(ft)) return false;
      if (fb && !(x.b || []).includes(fb)) return false;
      if (fd && !(x.d || []).some(v => v.toLowerCase().includes(fd))) return false;
      const fp = document.getElementById("cpath").value.trim().toLowerCase(),
            fo = document.getElementById("corigin").value,
            fr = document.getElementById("crole").value;
      if (fp && !(x.p || []).some(v => v.toLowerCase().includes(fp))) return false;
      if (fo && !(x.s || []).includes(fo)) return false;
      if (fr && !Object.values(x.l || {}).some(a => a.includes(fr))) return false;
      return true;
    }).sort((a, b) => (recCount[b[0]] || 0) - (recCount[a[0]] || 0));
    cpage = 0; crender();
  }
  function crender() {
    const pages = Math.max(1, Math.ceil(cview.length / PAGE));
    if (cpage >= pages) cpage = pages - 1;
    document.getElementById("ccount").textContent =
      `${cview.length} chemical${cview.length === 1 ? "" : "s"} of ${Object.keys(C).length}`;
    document.getElementById("crows").innerHTML =
      cview.slice(cpage * PAGE, cpage * PAGE + PAGE).map(([id, c]) => `<tr>
        <td class="pair">${esc(c.n)}${c.nt ? ' <span class="chip">neurotransmitter</span>' : ""}${c.pv ? ' <span class="chip">provisional</span>' : ""}</td>
        <td>${c.h ? `<a class="mono" href="https://hmdb.ca/metabolites/${esc(c.h)}" target="_blank" rel="noopener">${esc(c.h)}</a>` : '<span class="unassessed" title="name-keyed only">none</span>'}</td>
        <td class="mono" style="font-size:12px">${esc(c.f) || "&mdash;"}</td>
        <td>${c.mc ? `<span class="chip">${esc(c.mc)}</span>` : ""}</td>
        <td>${recCount[id] || 0}</td>
        <td>${new Set((E[id] || []).map(x => x[0])).size || '<span class="unassessed" title="not inferable from expression">0</span>'}</td>
        <td>${c.src.length}</td></tr>`).join("");
    document.getElementById("cpageinfo").textContent = `page ${cpage + 1} / ${pages}`;
    document.getElementById("cprev").disabled = cpage === 0;
    document.getElementById("cnext").disabled = cpage >= pages - 1;
  }

  /* ---------------- statistics ---------------- */
  function bars(obj, total) {
    const max = Math.max(...Object.values(obj), 1);
    return Object.entries(obj).map(([k, v]) => `<div class="hbar">
      <span class="lab">${esc(k)}</span>
      <span class="track"><i style="width:${(v / max * 100).toFixed(1)}%"></i></span>
      <span class="n">${v.toLocaleString()}</span></div>`).join("");
  }
  const TRUST = [
    ["GtoPdb (IUPHAR)", "Expert-curated pharmacology of receptors, with measured affinities", "Yes, per pair", "High: curated by receptor specialists; limited to well-studied drug targets"],
    ["CellChatDB v2", "Curated ligand-receptor list for cell-cell communication", "Some", "Medium-high: manually curated, partly built from KEGG; pair-level evidence is thin"],
    ["CellPhoneDB v5", "Curated ligand-receptor list, includes some non-protein ligands", "Rarely", "Medium-high: curated, but few citations per pair"],
    ["MRCLinkdb", "Metabolite-receptor pairs collected from the literature", "Yes, per pair", "Medium: citations present, but some are off by a few ids and the file has no stated licence"],
    ["MetaLigand", "Metabolite-ligand-receptor pairs, compiled from several databases", "No", "Low-medium: aggregation without citations; includes artefact-like pairs (for example glycerol with integrins)"],
    ["MetalinksDB", "Metabolite-protein links, mostly from other resources", "No", "Low-medium: about 94% of its provenance is STITCH (text-mining/prediction), so treat single-source pairs as leads"],
    ["ChEMBL", "Measured activities from medicinal-chemistry papers", "Yes, PubMed id", "High for 'binds/acts at' (Ki, Kd, EC50, IC50 up to 10 uM); species and assay vary"],
    ["BindingDB", "Measured binding affinities, matched by chemical structure", "Yes, PubMed id", "High for affinity; matched on the structure skeleton, so close analogues can slip in"],
    ["PDB", "3D structures that contain the molecule and a chain of the receptor", "Yes, primary paper", "Medium: shows co-occurrence in a structure, not that the molecule is the physiological ligand; common additives are excluded"],
    ["UniProt / Rhea", "Curated protein annotations (binding sites, reactions, function text)", "Yes, experimental PubMed id", "High when the sentence is experimental; reactions describe transport or catalysis as well as binding"],
    ["TCDB", "Transporter classification with substrates", "No", "Medium-high for transporters only"],
    ["HGNC names", "Receptor named after the ligand (for example 'adenosine A1 receptor')", "No", "High: the name itself says it"]
  ];
  function buildStats() {
    const m = S.manifest, single = S.by_nsrc["1"] || 0;
    const bs = S.by_source || {}, SRCL = { gtopdb_interactions: "GtoPdb", cellchatdb_human: "CellChatDB v2",
      mrclinkdb_human_mr: "MRCLinkdb", cellphonedb_interactions: "CellPhoneDB v5", metaligand_pairs_human: "MetaLigand",
      metalinksdb: "MetalinksDB" };
    const cmpMax = Math.max(S.n_inter, ...Object.values(bs).map(x => x.pairs), 1);
    const cmpBar = (lab, tot, ver, strong) => `<div class="hbar" title="${esc(lab)}: ${ver.toLocaleString()} of ${tot.toLocaleString()} pairs verified">
      <span class="lab"${strong ? ' style="color:var(--ink);font-weight:600"' : ""}>${esc(lab)}</span>
      <span class="track" style="position:relative"><i style="width:${(tot / cmpMax * 100).toFixed(1)}%;opacity:.35"></i>
        <i style="position:absolute;left:0;top:0;width:${(ver / cmpMax * 100).toFixed(1)}%;opacity:.9"></i></span>
      <span class="n" style="width:150px">${tot.toLocaleString()} &middot; ${ver.toLocaleString()} verified</span></div>`;
    const cmpRows = Object.entries(bs).map(([k, v]) => cmpBar(SRCL[k] || k, v.pairs, v.own_verified ?? 0, false)).join("") +
      cmpBar("This database", S.n_inter, S.n_verified || 0, true) +
      `<div class="donut-key"><span><i style="background:var(--accent)"></i>verified</span>
        <span><i style="background:var(--accent);opacity:.35"></i>listed, not verified</span></div>`;
    document.getElementById("statbody").innerHTML = `
      <div class="stat-grid">
        <div class="stat"><b>${S.n_inter.toLocaleString()}</b><span>chemical–receptor pairs</span></div>
        <div class="stat"><b>${S.n_chem}</b><span>chemical ligands</span></div>
        <div class="stat"><b>${(S.n_verified || 0).toLocaleString()}</b><span>verified pairs</span></div>
      </div>
      <div class="panel" style="margin-bottom:16px">
        <h3 style="margin:0 0 8px;font-size:15px">Existing resources compared with this database</h3>
        <p style="color:var(--muted);font-size:12.5px;margin:0 0 10px">Each existing resource as it stands: how many
          chemical&ndash;receptor pairs it lists, and how many of those it already verifies by this project's standard, using only
          what the resource itself supplies (a paper it cites that states the interaction in a quoted sentence, or a curated
          record with measured pharmacology, as in GtoPdb). The last bar is this database, after adding literature reading and
          curated-database records (ChEMBL, BindingDB, PDB, UniProt/Rhea, TCDB) across all pairs.</p>
        ${cmpRows}
      </div>
      <div class="panel" style="margin-bottom:16px">
        <h3 style="margin:0 0 8px;font-size:15px">Database evidence behind verified pairs</h3>
        <p style="color:var(--muted);font-size:12.5px;margin:0 0 10px">${(S.n_lit_verified || 0).toLocaleString()} pairs
          are verified by a quoted sentence from a paper; the database records below verify further pairs
          (a pair can appear under several). These are measured activities or structures recorded by the
          database, not quoted sentences.</p>
        ${bars(S.by_dbsource || {})}
      </div>
      <div class="panel" style="margin-bottom:16px">
        <h3 style="margin:0 0 8px;font-size:15px">Enzyme roles</h3>
        <p style="color:var(--muted);font-size:12.5px;margin:0 0 10px">The five-role union
          vocabulary, mapped from five different upstream schemes.</p>
        ${bars(S.by_role)}
      </div>
      <div class="panel" style="margin-bottom:16px">
        <h3 style="margin:0 0 8px;font-size:15px">Molecular class </h3>
        ${bars(S.by_class)}
      </div>
      <div class="panel">
        <h3 style="margin:0 0 8px;font-size:15px">Receptor type</h3>
        ${bars(S.by_rtype)}
      </div>`;
  }

  /* ---------------- methods ---------------- */
  function buildMethods() {
    const m = S.manifest;
    const tbl = (head, rows) => `<table class="trust"><thead><tr>${head.map(h => `<th>${h}</th>`).join("")}</tr></thead><tbody>${
      rows.map(r => `<tr>${r.map(x => `<td>${x}</td>`).join("")}</tr>`).join("")}</tbody></table>`;
    document.getElementById("methodsbody").innerHTML = `
      <h3 style="margin-top:0">1. Databases the pairs were aggregated from</h3>
      <p>Every chemical&ndash;receptor pair in the table comes from at least one of these. They decide
      <i>which pairs exist</i> here; listing a pair does not make it verified.</p>
      ${tbl(["Database", "What it contributes", "Pairs"], [
        ["GtoPdb (IUPHAR/BPS Guide to Pharmacology)", "Expert-curated receptor pharmacology with affinities and citations", (S.by_source || {}).gtopdb_interactions?.pairs ?? ""],
        ["CellChatDB v2", "Curated ligand&ndash;receptor list for cell&ndash;cell communication (non-protein ligands)", (S.by_source || {}).cellchatdb_human?.pairs ?? ""],
        ["CellPhoneDB v5", "Curated ligand&ndash;receptor list, non-protein ligands only here", (S.by_source || {}).cellphonedb_interactions?.pairs ?? ""],
        ["MRCLinkdb", "Metabolite&ndash;receptor pairs from the literature, with PubMed ids", (S.by_source || {}).mrclinkdb_human_mr?.pairs ?? ""],
        ["MetaLigand", "Metabolite ligand&ndash;receptor pairs compiled from other databases", (S.by_source || {}).metaligand_pairs_human?.pairs ?? ""],
        ["MetalinksDB", "Metabolite&ndash;protein links, largely from STITCH", (S.by_source || {}).metalinksdb?.pairs ?? ""],
        ["MEBOCOST", "Metabolite&ndash;sensor pairs; used internally only, never shown or exported (licence)", "&mdash;"]])}
      <p style="color:var(--muted);font-size:12.5px">The genes for each ligand (synthesis, precursor synthesis, export,
      uptake, breakdown) are taken from the raw sources rather than from the ligand&ndash;receptor resources: the
      Human-GEM metabolic model (v2.1.1) and Rhea reactions as annotated in UniProt, using only reactions with a stated
      direction. For the few chemicals these sources do not cover, the genes listed by MetalinksDB, MetaLigand or MRCLinkdb
      are shown instead and marked as not confirmed. Cofactors such as ATP, ADP, phosphate and NAD get no gene set: hundreds
      of enzymes make or use them. Gene symbols are current HGNC names.</p>

      <h3>2. Databases and literature used to verify pairs</h3>
      <p>A pair is <b>verified</b> when one of these backs it independently of the list it came from.</p>
      ${tbl(["Source", "What counts as evidence", "Pairs it verifies"], [
        ["Primary papers (PubMed, Europe PMC, PMC full text, publisher pages)", "A sentence from the paper's own Results or Discussion, quoted word for word, naming this molecule and this receptor", (S.n_lit_verified || 0).toLocaleString()],
        ...Object.entries(S.by_dbsource || {}).map(([k, v]) => [esc(k), ({
          ChEMBL: "A measured activity (Ki, Kd, IC50, EC50 at 10 &micro;M or stronger) on the human protein, with a PubMed id",
          BindingDB: "A measured binding affinity for the same chemical structure on this protein",
          GtoPdb: "An expert-curated interaction record (agonist, antagonist, affinity)",
          PDB: "A 3D structure containing the molecule and a chain of the receptor (common additives excluded)",
          UniProt: "A curated binding-site or function annotation with an experimental PubMed id",
          "UniProt (Rhea)": "A curated reaction (transport or catalysis) involving the molecule",
          TCDB: "A transporter substrate annotation",
          "HGNC approved name": "The receptor is officially named after the molecule"
        })[k] || "", v.toLocaleString()])])}
      <p style="color:var(--muted);font-size:12.5px">A pair can be verified by several of these at once.</p>`;
  }


  /* ---------------- wiring ---------------- */
  function initSelects() {
    const cls = [...new Set(Object.values(C).map(c => c.mc).filter(Boolean))].sort();
    document.getElementById("mclass").insertAdjacentHTML("beforeend",
      cls.map(c => `<option>${esc(c)}</option>`).join(""));
    /* Assay classes come from the confirmed evidence rows. Each option shows how many
       pairs in the table have at least one evidence row of that class. */
    const ASSAY_LABEL = { structural: "structural", direct_binding: "direct binding",
      functional_in_vitro: "functional in vitro", functional_in_vivo: "functional in vivo",
      computational: "computational" };
    const perClass = {};
    I.forEach(r => {
      new Set(((window.EVIDENCE || {})[r[0]] || []).map(e => e.as).filter(Boolean))
        .forEach(k => { perClass[k] = (perClass[k] || 0) + 1; });
    });
    document.getElementById("assay").insertAdjacentHTML("beforeend",
      Object.keys(ASSAY_LABEL).filter(k => perClass[k]).map(k =>
        `<option value="${k}">${esc(ASSAY_LABEL[k])} (${perClass[k]})</option>`).join(""));
    const tally = k => { const t = {}; Object.values(C).forEach(c => ((c.x || {})[k] || []).forEach(v => { t[v] = (t[v] || 0) + 1; }));
      return Object.entries(t).sort((p, q) => q[1] - p[1] || p[0].localeCompare(q[0])); };
    const tallyRole = () => { const t = {}; Object.values(C).forEach(c => Object.values((c.x || {}).l || {})
      .forEach(a => new Set(a).forEach(v => { t[v] = (t[v] || 0) + 1; })));
      return Object.entries(t).sort((p, q) => q[1] - p[1] || p[0].localeCompare(q[0])); };
    document.getElementById("corigin").insertAdjacentHTML("beforeend",
      tally("s").map(([v, n]) => `<option value="${esc(v)}">${esc(v)} (${n})</option>`).join(""));
    document.getElementById("crole").insertAdjacentHTML("beforeend",
      tallyRole().map(([v, n]) => `<option value="${esc(v)}">${esc(v)} (${n})</option>`).join(""));
    [["ctissue", "t"], ["cbio", "b"]].forEach(([id, k]) =>
      document.getElementById(id).insertAdjacentHTML("beforeend",
        tally(k).map(([v, n]) => `<option value="${esc(v)}">${esc(v)} (${n})</option>`).join("")));
    const rts = [...new Set(I.map(r => r[4]).filter(Boolean))].sort();
    document.getElementById("rtype").insertAdjacentHTML("beforeend",
      rts.map(c => `<option>${esc(c)}</option>`).join(""));

  }

  document.querySelectorAll("nav button").forEach(b => b.onclick = () => {
    document.querySelectorAll("nav button").forEach(x => x.classList.remove("on"));
    b.classList.add("on");
    ["browse", "chems", "stats", "methods", "usage"].forEach(t =>
      document.getElementById("tab-" + t).classList.toggle("hide", t !== b.dataset.tab));
    if (b.dataset.tab === "chems" && !cview.length) capply();
  });
  ["q", "mclass", "rtype", "minsrc", "assay", "only"].forEach(id => {
    const el = document.getElementById(id);
    el.addEventListener(id === "q" ? "input" : "change", apply);
  });
  ["cq", "conly", "ctissue", "cbio", "cdis", "cpath", "corigin", "crole"].forEach(id => document.getElementById(id)
    .addEventListener(id === "cq" || id === "cdis" || id === "cpath" ? "input" : "change", capply));
  document.querySelectorAll("th[data-s]").forEach(th => th.onclick = () => {
    const c = +th.dataset.s;
    if (c === sortCol) sortAsc = !sortAsc; else { sortCol = c; sortAsc = false; }
    document.querySelectorAll("th[data-s]").forEach(x => x.classList.remove("sorted", "asc"));
    th.classList.add("sorted"); if (sortAsc) th.classList.add("asc");
    apply();
  });
  document.getElementById("dlcsv").onclick = downloadChemCSV;
  document.getElementById("prev").onclick = () => { page--; render(); };
  document.getElementById("next").onclick = () => { page++; render(); };
  document.getElementById("cprev").onclick = () => { cpage--; crender(); };
  document.getElementById("cnext").onclick = () => { cpage++; crender(); };
  document.getElementById("rows").onclick = e => {
    const tr = e.target.closest("tr[data-id]");
    if (tr) openDetail(tr.dataset.id);
  };
  document.addEventListener("keydown", e => {
    if (e.key === "Escape") document.getElementById("detail").innerHTML = "";
  });

  initSelects(); apply(); buildStats(); buildMethods();
})();
