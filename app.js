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
    metaligand_ligands_human: "MetaLigand", mrclinkdb_human_enzyme: "MRCLinkdb", mrclinkdb_human_transporter: "MRCLinkdb"
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
  function enzPanel(cid) {
    const rows = E[cid] || [], GI = window.GENEINFO || {};
    if (!rows.length) return `<p style="color:var(--muted)">No genes recorded. This chemical
      <b>cannot be inferred from transcriptomic data</b> — there is no gene whose expression
      stands in for it.</p>`;
    const byRole = {};
    rows.forEach(([g, role, src, dk]) => {
      (byRole[role] = byRole[role] || {});
      (byRole[role][g] = byRole[role][g] || { srcs: new Set(), dk: dk });
      byRole[role][g].srcs.add(src);
    });
    const order = ["precursor_synthesis", "synthesis", "efflux", "influx", "degradation"];
    return `<p style="color:var(--muted);font-size:12px;margin:0 0 4px">Within each role, genes are ordered by how many
      source databases list them for this chemical (the number on each gene); the tag says what kind of protein it is.</p>` +
      order.filter(r => byRole[r]).map(role => {
      const genes = Object.entries(byRole[role]).sort((a, b) => b[1].srcs.size - a[1].srcs.size || (a[0] < b[0] ? -1 : 1));
      return `<div style="margin:8px 0">
        <div style="font-size:12px;color:var(--muted)" title="${esc(ROLE_HELP[role] || "")}">
          ${esc(role.replace("_", " "))} · ${genes.length} gene${genes.length === 1 ? "" : "s"}</div>
        <div style="margin-top:4px">${genes.map(([g, v]) => {
          const gi = GI[g] || ["other", ""], t = gi[0];
          return `<span class="pill ${ROLE_CLASS[role] || ""}" title="${esc((gi[1] ? gi[1] + " — " : "") + "listed by " + [...v.srcs].map(s => SRC_LABEL[s] || s).join(", "))}${v.dk ? "" : " — direction not stated by the source"}">${esc(g)}${v.dk ? "" : "?"}<sup class="gcount">${v.srcs.size}</sup><span class="gtype gt-${t}" title="${esc(GT_HELP[t])}">${esc(GT_LABEL[t])}</span></span>`;
        }).join("")}</div></div>`;
    }).join("");
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
          ${ch.sm ? `<dt title="SMILES is a one-line text code for a molecule's structure: letters are atoms, brackets and digits show branches and rings. Software can redraw the structure from it.">SMILES <span style="color:var(--muted)" title="A one-line text code for the molecule's structure (atoms and bonds); software can redraw the structure from it.">(?)</span></dt><dd class="smiles">${esc(ch.sm)}</dd>` : ""}
          ${ch.syn && ch.syn.length ? `<dt>Also known as</dt><dd style="font-size:12.5px;color:var(--muted)">${esc(ch.syn.slice(0, 8).join(", "))}</dd>` : ""}
        </dl>

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
    const cmpBar = (lab, tot, ver, strong) => `<div class="hbar">
      <span class="lab"${strong ? ' style="color:var(--text,#111);font-weight:600"' : ""}>${esc(lab)}</span>
      <span class="track" style="position:relative"><i style="width:${(tot / cmpMax * 100).toFixed(1)}%;opacity:.35"></i>
        <i style="position:absolute;left:0;top:0;width:${(ver / cmpMax * 100).toFixed(1)}%;opacity:.9"></i></span>
      <span class="n" style="width:120px">${tot.toLocaleString()} &middot; ${ver.toLocaleString()} verified</span></div>`;
    const cmpRows = Object.entries(bs).map(([k, v]) => cmpBar(SRCL[k] || k, v.pairs, v.verified, false)).join("") +
      cmpBar("This database", S.n_inter, S.n_verified || 0, true);
    document.getElementById("statbody").innerHTML = `
      <div class="stat-grid">
        <div class="stat"><b>${S.n_inter.toLocaleString()}</b><span>chemical–receptor pairs</span></div>
        <div class="stat"><b>${S.n_chem}</b><span>chemical ligands</span></div>
        <div class="stat"><b>${(S.n_verified || 0).toLocaleString()}</b><span>verified pairs</span></div>
      </div>
      <div class="panel" style="margin-bottom:16px">
        <h3 style="margin:0 0 8px;font-size:15px">Pairs in each resource, and how many we verified</h3>
        <p style="color:var(--muted);font-size:12.5px;margin:0 0 10px">Number of chemical&ndash;receptor pairs each
          upstream resource lists (light bar), and how many of those pairs this database could verify
          (dark part). The last row is this database: every pair it holds, and the pairs it
          verified (literature with a checked quote, or a curated-database record).</p>
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
        <h3 style="margin:0 0 8px;font-size:15px">Where the data comes from, and how far to trust it</h3>
        <table class="trust"><thead><tr><th>Source</th><th>What it is</th><th>Cites papers?</th><th>How much to trust it</th></tr></thead><tbody>
          ${TRUST.map(t => `<tr><td><b>${esc(t[0])}</b></td><td>${esc(t[1])}</td><td>${esc(t[2])}</td><td>${esc(t[3])}</td></tr>`).join("")}
        </tbody></table>
        <p style="color:var(--muted);font-size:12px;margin:10px 0 0">Trust notes are this project's own assessment.
          A pair listed by only one text-mined or prediction-based resource is a lead, not a fact; it
          becomes <i>verified</i> here only when a paper's own sentence or a curated database record backs it.</p>
      </div>
      <div class="panel" style="margin-bottom:16px">
        <h3 style="margin:0 0 8px;font-size:15px">How many sources agree</h3>
        <p style="color:var(--muted);font-size:12.5px;margin:0 0 10px">
          ${single.toLocaleString()} of ${S.n_inter.toLocaleString()} pairs
          (${(single / S.n_inter * 100).toFixed(0)}%) are carried by a single source; the rest
          are reported by two or more.</p>
        ${bars(S.by_nsrc)}
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
      <p style="color:var(--muted);font-size:12.5px">The genes for each ligand (synthesis, transport, degradation)
      come from MetalinksDB, MetaLigand, MRCLinkdb, CellChatDB and CellPhoneDB, with gene symbols brought to the
      current HGNC names.</p>

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
      <p style="color:var(--muted);font-size:12.5px">A pair can be verified by several of these at once.</p>

      <h3>3. What "verified" means here</h3>
      <ul>
        <li><b>Verified</b>: at least one quoted sentence from a primary paper, or at least one curated database record from section 2.</li>
        <li><b>Not verified</b>: only the aggregated lists in section 1 assert the pair. Treat it as a lead.</li>
        <li>A paper that was read but does not test the pair is <i>not</i> counted as evidence against it. Only a paper that tested
          this exact pair and found nothing counts as a refutation.</li>
        <li>Quoted sentences were checked mechanically against the paper's own text; they were not read by a person.</li>
      </ul>

      <h3>What this is</h3>
      <p>Seven existing databases describe interactions between small molecules —
      metabolites, neurotransmitters, steroids, lipids — and the proteins that receive them.
      Each uses its own names, its own identifiers and its own table layout. This resource
      merges all seven into one table, one row per chemical–receptor pair in human, and keeps
      every original row alongside it so any entry can be traced back to the database it
      came from.</p>
      <p>It exists because the peptide resource is built on connectomeDB2025, which excludes
      non-peptide ligands as a matter of policy. Glutamate, GABA, serotonin, acetylcholine,
      the prostaglandins and the steroids are absent from it entirely.</p>

      <h3>The sources</h3>
      <div>${(m.sources || []).filter(s => s.shippable).map(s =>
        `<div class="dl"><b>${esc(s.name)}</b><span>${esc(s.licence)}</span></div>`).join("")}</div>

      <h3>How they were merged, and the choices that involved</h3>

      <p><b>Matching the chemicals.</b> The hard part is knowing when two databases mean the
      same molecule. They rarely spell it the same way: one writes
      <code>7-a,25-Dihydroxycholesterol</code>, another <code>7α,25-dihydroxycholesterol</code>.
      Four of the seven give an HMDB identifier, so those were matched on the identifier.
      The rest were matched by name, then by a list of synonyms, in that order. Where two
      sources disagreed about a molecule's identity the disagreement was recorded rather than
      resolved silently, and the first value seen was kept.</p>

      <p><b>Chemicals were given their own table.</b> A metabolite has no gene and no protein
      accession, so it cannot sit in a table built for genes without inventing something. It
      gets its own record with HMDB and PubChem identifiers, a formula and a structure where
      a source supplied them. ChEBI identifiers are empty throughout: no source in the union
      publishes one, and none was guessed.</p>

      <p><b>Two sources bundle the ligand and its enzymes into one label.</b> CellChatDB
      writes <code>5-HT-TPH1_SLC18A1</code>, meaning serotonin, made by TPH1, exported by
      SLC18A1. That is three separate facts in one string. All 239 such labels were taken
      apart so the chemical is stored as the chemical and the enzymes are stored separately.
      Fifteen needed a person to decide; the rest were unambiguous.</p>

      <p><b>Enzymes are kept apart from receptors, deliberately.</b> You cannot measure a
      metabolite directly in transcriptomic data — it has no transcript. What you can measure
      are the genes that make it, export it, import it or break it down. Those genes are
      recorded in a separate table, because "this cell can produce the molecule" is a
      different statement from "this molecule binds that receptor", and merging the two would
      blur them. Five roles are used: synthesis, precursor synthesis, degradation, export and
      import. Each source labels these differently, so every label was mapped onto the five
      and the original wording kept next to it.</p>

      <p><b>Where sources disagree, both readings are kept.</b> For dopamine, three databases
      call SLC18A1 and SLC18A2 exporters and one calls them importers. Both are defensible —
      the transporter moves dopamine into a vesicle, which is import across that membrane and
      export from the cell's point of view. Rather than pick, the table carries both with the
      source attached.</p>

      <p><b>Receptor complexes were broken into their parts.</b> Three of the sources describe
      multi-protein receptors, each in its own notation — <code>GABR_A1B2</code>,
      <code>GABRA1;GABRG2;GABRB2</code>, and a named entry in a separate file. Each was split
      into individual proteins, with a label tying the parts back together.</p>

      <p><b>Some rows were dropped.</b> Eight entries filed under non-peptide signalling are
      peptides or protein complexes, and were excluded. A few rows fused two different
      chemicals into one entry — one source lists "Aldosterone Corticosterone" as a single
      ligand, another pairs glycerol with fumaric acid — and those were split apart by
      identifier rather than by splitting the text, because several chemical names contain
      commas themselves.</p>

      <p><b>One source is not redistributed.</b> MEBOCOST's licence does not permit it, so its
      rows are withheld and its name is removed from the source lists shown here. Counts of
      how many sources report a pair can therefore be lower here than in the working
      database.</p>

      <h3>Confirmed evidence, and what that word means here</h3>
      <p>${S.n_evidence || 0} primary papers have been read for ${S.n_pairs_with_evidence || 0}
      of the ${S.n_inter.toLocaleString()} pairs, each one quoting the sentence in the paper that
      states the result. Those pairs carry a green panel in the drawer, and the
      <i>has confirmed evidence</i> filter shows only them. Everything else carries provenance
      alone.</p>
      <p><b>Confirmed here means machine-verified, not curator-read.</b> Each quoted sentence was
      checked to appear verbatim in the paper's own abstract or publisher full text, with the
      molecule and the receptor both named in that same text, and the paper checked not to be a
      review or retracted. No person has read the paper. That is a weaker claim than the peptide
      half's confirmations, and the two should not be compared.</p>
      <p>The labels answer different questions. <b>Supports</b> and <b>refutes</b> are about the
      pair as stated: a refutation means a paper tested this molecule at this receptor and found
      nothing, which makes the pair wrong rather than merely unevidenced. There are
      ${(S.by_result || {}).REFUTES || 0} of them. <b>Family level</b> means the paper implicates
      the receptor family but not the named subunit. <b>Transporter</b> and <b>enzyme</b> mean the
      protein moves or consumes the molecule, which is a real relationship but not receptor
      binding.</p>
      <p>Absence of a panel is not evidence of absence. Of 435 pairs searched so far, 298 returned
      no qualifying paper, and the reasons differ sharply: the wrong molecule (cystine where the
      pair says cysteine, cholecalciferol where the receptor wants calcitriol), the cargo mistaken
      for the ligand (iron at TFR1, whose ligand is transferrin), a dimerisation partner mistaken
      for the target (9-cis-retinoic acid at every RXR partner), a name collision (L-arginine at
      the arginine-vasopressin receptor), a human pseudogene (GLRA4), or simply nobody having run
      the experiment. Those reasons are recorded per pair in the working database, not here.</p>

      <h3>Reading a row</h3>
      <p>Each row names a chemical and a receptor, the sources that report the pair, and the
      PubMed identifiers those sources cite, where they cite any. Clicking a row opens the
      chemical's identifiers and classification, the genes that make and move it grouped by
      role, any confirmed evidence, and every original database row behind the pair, exactly as
      that database wrote it.</p>`;
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
  ["cq", "conly"].forEach(id => document.getElementById(id)
    .addEventListener(id === "cq" ? "input" : "change", capply));
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
