/* Usage tab: how the five enzyme/transporter roles can be turned into an
   estimate of how available a ligand is in a cell type or sample.

   The formulation is MetaLigand's (Cell Rep Methods 2025, PMID 41175863): de novo and
   salvage synthesis, efflux and receptor, combined with square roots of means (the genes of
   one role are averaged, the sequential steps are multiplied under a square root).
   Degradation genes are NOT used: subtracting them (an earlier version of this tab, after
   MRCLinkdb's producers-minus-consumers idea) lowered the ranking of known producer tissues
   and cell types in our benchmark (formula_benchmark/RESULTS.md). The formula is not
   validated against measured concentrations. */
(function () {
  "use strict";
  const esc = s => String(s == null ? "" : s).replace(/[&<>"']/g,
    m => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[m]));

  const mean = a => a.length ? a.reduce((x, y) => x + y, 0) / a.length : null;
  const f3 = x => (x == null ? "n/a" : (Math.round(x * 1000) / 1000).toString());

  /* The formulation. Inputs are arrays of non-negative normalised expression values,
     one value per gene in that role. Returns every intermediate so the page can show it. */
  function score(v) {
    const syn = mean(v.synthesis), pre = mean(v.precursor), inf = mean(v.influx),
          eff = mean(v.efflux), rec = mean(v.receptor);
    const r = { syn, pre, inf, eff, rec, steps: [] };
    if (syn == null) { r.error = "Needs at least one synthesis gene."; return r; }
    r.SD = pre != null ? Math.sqrt(syn * pre) : syn;
    r.steps.push(pre != null
      ? `S_D = sqrt(mean(synthesis) x mean(precursor)) = sqrt(${f3(syn)} x ${f3(pre)}) = ${f3(r.SD)}`
      : `S_D = mean(synthesis) = ${f3(syn)}   (no precursor genes entered)`);
    r.SS = inf != null ? Math.sqrt(syn * inf) : null;
    if (r.SS != null) r.steps.push(`S_S = sqrt(mean(synthesis) x mean(influx)) = sqrt(${f3(syn)} x ${f3(inf)}) = ${f3(r.SS)}`);
    r.S = r.SS != null ? (r.SD + r.SS) / 2 : r.SD;
    r.steps.push(r.SS != null
      ? `S = (S_D + S_S) / 2 = ${f3(r.S)}`
      : `S = S_D = ${f3(r.S)}   (no influx genes entered, so no salvage term)`);
    r.M = eff != null ? Math.sqrt(r.S * eff) : r.S;
    r.steps.push(eff != null
      ? `M = sqrt(S x mean(efflux)) = sqrt(${f3(r.S)} x ${f3(eff)}) = ${f3(r.M)}`
      : `M = S = ${f3(r.M)}   (no efflux genes entered)`);
    if (rec != null) {
      r.L = Math.sqrt(r.M * rec);
      r.steps.push(`Signal = sqrt(M x mean(receptor)) = sqrt(${f3(r.M)} x ${f3(rec)}) = ${f3(r.L)}`);
    }
    return r;
  }

  const EXAMPLE = {
    names: {
      synthesis: "DDC", precursor: "TH", influx: "(none in the core set)",
      efflux: "SLC18A2, SLC6A3", receptor: "DRD2"
    },
    cells: {
      "Dopaminergic neuron (illustrative)": {
        synthesis: [4.0], precursor: [3.6], influx: [],
        efflux: [3.0, 2.4], receptor: [2.5]
      },
      "Non-dopaminergic cell (illustrative)": {
        synthesis: [3.0], precursor: [0.1], influx: [],
        efflux: [0.2, 0.0], receptor: [2.5]
      }
    }
  };

  function stepsHtml(r) {
    return r.error ? `<p style="color:#8a2b2b">${esc(r.error)}</p>` :
      `<pre class="mono" style="background:var(--chip);padding:10px;border-radius:6px;overflow-x:auto;font-size:12.5px;margin:6px 0">${esc(r.steps.join("\n"))}</pre>`;
  }

  function build() {
    const first = Object.values(EXAMPLE.cells)[0];
    const box = (k, label) => `<div style="margin:6px 0"><label style="font-size:12.5px;display:block">${label}</label>
      <input type="text" id="uc-${k}" value="${esc(first[k].join(", "))}" style="width:100%;max-width:420px"></div>`;

    document.getElementById("usagebody").innerHTML = `
      <h3 style="margin-top:0">Estimating a ligand from enzyme and transporter genes</h3>
      <p>A metabolite, steroid or lipid has no transcript, so it cannot be seen directly in
      single-cell or bulk RNA data. What can be measured is the expression of the genes that make,
      remove and move it. This tab shows how the five enzyme roles on this site can be combined
      into one number for a cell type or sample, and how that number feeds a ligand&ndash;receptor
      score. It is an <b>estimate of availability from RNA</b>, not a measured concentration.</p>

      <h3>The gene roles</h3>
      <ul style="font-size:13.5px;line-height:1.7">
        <li><b>Precursor synthesis</b>: genes that make the molecule the ligand is built from.</li>
        <li><b>Synthesis</b>: genes that make the ligand itself.</li>
        <li><b>Influx</b>: transporters that bring the ligand, or its precursor, into the cell to be reused.</li>
        <li><b>Efflux</b>: transporters that move the ligand out of the cell or load it into vesicles for release.</li>
        <li><b>Degradation</b>: genes that break the ligand down or inactivate it. Listed on this site,
        but <b>not used</b> in the formula below (see the note after it).</li>
      </ul>

      <h3>Which genes to use</h3>
      <p style="font-size:13.5px">The gene lists on this site come from the raw sources (the Human-GEM metabolic model
      and Rhea/UniProt reactions, directed reactions only), so the roles below can be used as given. Faded genes, shown
      only for chemicals the raw sources do not cover, are not confirmed: use them with care. <b>Cofactor hubs</b>
      (ATP, ADP, phosphate, NAD and similar) have no gene list: hundreds of enzymes make or use them, so no gene set
      can stand in for where they are released as signals.</p>

      <h3>The formula</h3>
      <pre class="mono" style="background:var(--chip);padding:12px;border-radius:6px;overflow-x:auto;font-size:12.5px;line-height:1.7">1  Production
     S_D = sqrt( mean(synthesis) x mean(precursor) )
     S_S = sqrt( mean(synthesis) x mean(influx) )
     S   = (S_D + S_S) / 2

2  Availability
     M   = sqrt( S x mean(efflux) )

3  Signal strength for one ligand-receptor pair
     Signal = sqrt( M x mean(receptor genes) )</pre>
      <p style="font-size:13.5px">This is the formula of MetaLigand (Cell Reports Methods 2025). The genes of one
      role are alternatives and are averaged; the steps that follow each other are multiplied under a square root,
      so a missing step pulls the score down. If a role has no genes, its step is skipped: without precursor genes
      S_D = mean(synthesis); without influx genes (no salvage route) S = S_D; without efflux genes M = S.
      At least one synthesis gene is needed.</p>
      <p style="font-size:13.5px"><b>Why degradation is not subtracted.</b> In a test of 15 formulas on the same gene
      sets, ranking the known producer tissues and cell types of 37 and 44 ligands (Human Protein Atlas), this
      formula scored a mean AUROC of 0.878 (tissues) and 0.858 (cell types), against 0.836 and 0.738 for the same
      formula with degradation subtracted. Subtracting sets a producer to zero when its breakdown genes are also
      high. The test shows where a ligand is made, not how much; it has not been checked against measured
      concentrations.</p>

      <h3>Try your own numbers</h3>
      <p style="font-size:12.5px;color:var(--muted)">Enter normalised expression values for each
      role, separated by commas (one value per gene). Leave a box empty if the ligand has no genes
      in that role.</p>
      ${box("synthesis", "Synthesis genes")}
      ${box("precursor", "Precursor synthesis genes")}
      ${box("influx", "Influx genes (salvage)")}
      ${box("efflux", "Efflux genes")}
      ${box("receptor", "Receptor genes")}
      <p><button id="uc-run">Calculate</button></p>
      <div id="uc-out"></div>`;

    const read = id => document.getElementById("uc-" + id).value
      .split(/[\s,;]+/).filter(Boolean).map(Number).filter(x => isFinite(x) && x >= 0);
    const KEYS = ["synthesis", "precursor", "influx", "efflux", "receptor"];
    const run = () => {
      const v = {}; KEYS.forEach(k => v[k] = read(k));
      document.getElementById("uc-out").innerHTML = stepsHtml(score(v));
    };
    document.getElementById("uc-run").onclick = run;
    run();
  }

  window.buildUsage = build;
  build();
})();
