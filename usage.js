/* Usage tab: how the five enzyme/transporter roles can be turned into an
   estimate of how available a ligand is in a cell type or sample.

   The formulation is OURS, adapted from two published ones:
     MetaLigand (Cell Rep Methods 2025, PMID 41175863): de novo / salvage synthesis,
       efflux, receptor, combined with square roots of means.
     MRCLinkdb (BMC Biol 2024, PMID 38978014): producers minus consumers.
   Neither paper publishes this exact combination; it is not validated. */
(function () {
  "use strict";
  const esc = s => String(s == null ? "" : s).replace(/[&<>"']/g,
    m => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[m]));
  const S = window.STATS || {};

  const mean = a => a.length ? a.reduce((x, y) => x + y, 0) / a.length : null;
  const f3 = x => (x == null ? "n/a" : (Math.round(x * 1000) / 1000).toString());

  /* The formulation. Inputs are arrays of non-negative normalised expression values,
     one value per gene in that role. Returns every intermediate so the page can show it. */
  function score(v) {
    const syn = mean(v.synthesis), pre = mean(v.precursor), inf = mean(v.influx),
          deg = mean(v.degradation), eff = mean(v.efflux), rec = mean(v.receptor);
    const r = { syn, pre, inf, deg, eff, rec, steps: [] };
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
    const d = deg == null ? 0 : deg;
    r.N = Math.max(0, r.S - d);
    r.steps.push(deg != null
      ? `N = max(0, S - mean(degradation)) = max(0, ${f3(r.S)} - ${f3(deg)}) = ${f3(r.N)}`
      : `N = S = ${f3(r.N)}   (no degradation genes entered)`);
    r.M = eff != null ? Math.sqrt(r.N * eff) : r.N;
    r.steps.push(eff != null
      ? `M = sqrt(N x mean(efflux)) = sqrt(${f3(r.N)} x ${f3(eff)}) = ${f3(r.M)}`
      : `M = N = ${f3(r.M)}   (no efflux genes entered)`);
    if (rec != null) {
      r.L = Math.sqrt(r.M * rec);
      r.steps.push(`Signal = sqrt(M x mean(receptor)) = sqrt(${f3(r.M)} x ${f3(rec)}) = ${f3(r.L)}`);
    }
    return r;
  }

  const EXAMPLE = {
    names: {
      synthesis: "DDC", precursor: "TH", influx: "(none in the core set)",
      degradation: "COMT, MAOA, MAOB", efflux: "SLC18A2, SLC6A3", receptor: "DRD2"
    },
    cells: {
      "Dopaminergic neuron (illustrative)": {
        synthesis: [4.0], precursor: [3.6], influx: [], degradation: [1.5, 0.8, 1.0],
        efflux: [3.0, 2.4], receptor: [2.5]
      },
      "Non-dopaminergic cell (illustrative)": {
        synthesis: [3.0], precursor: [0.1], influx: [], degradation: [3.0, 2.5, 3.5],
        efflux: [0.2, 0.0], receptor: [2.5]
      }
    }
  };

  const ROLES = [
    ["synthesis", "Synthesis", "Genes whose enzymes make the ligand itself",
     "MetalinksDB, MetaLigand, MRCLinkdb (producers), CellChatDB, CellPhoneDB"],
    ["precursor_synthesis", "Precursor synthesis", "Genes that make the ligand's precursor, not the ligand",
     "MetaLigand only"],
    ["influx", "Influx", "Transporters that bring the ligand or its precursor into the cell (salvage)",
     "MetalinksDB, MetaLigand"],
    ["degradation", "Degradation", "Genes whose enzymes consume or inactivate the ligand",
     "MetalinksDB, MRCLinkdb (consumers)"],
    ["efflux", "Efflux", "Transporters that move the ligand out of the cell, including vesicle loading",
     "MetalinksDB, MetaLigand, MRCLinkdb (no direction stated), CellChatDB, CellPhoneDB"]
  ];

  function stepsHtml(r) {
    return r.error ? `<p style="color:#8a2b2b">${esc(r.error)}</p>` :
      `<pre class="mono" style="background:var(--chip);padding:10px;border-radius:6px;overflow-x:auto;font-size:12.5px;margin:6px 0">${esc(r.steps.join("\n"))}</pre>`;
  }

  function build() {
    const by = S.by_role || {};
    const roleRows = ROLES.map(([k, label, what, src]) =>
      `<tr><td><b>${esc(label)}</b></td><td>${esc(what)}</td><td>${esc(src)}</td>
       <td style="text-align:right">${(by[k] || 0).toLocaleString()}</td></tr>`).join("");

    const exCells = Object.entries(EXAMPLE.cells).map(([name, v]) => {
      const r = score(v);
      return `<h4 style="margin:14px 0 4px;font-size:13.5px">${esc(name)}</h4>
        <div style="font-size:12.5px;color:var(--muted)">expression:
          ${["synthesis","precursor","degradation","efflux","receptor"].map(k =>
            `${k} ${esc(EXAMPLE.names[k])} = ${esc(v[k].join(", "))}`).join(" &middot; ")}</div>
        ${stepsHtml(r)}`;
    }).join("");

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

      <h3>The five roles and where they come from</h3>
      <div style="overflow-x:auto"><table>
        <thead><tr><th>Role</th><th>What the genes do</th><th>Source databases</th><th style="text-align:right">Links on this site</th></tr></thead>
        <tbody>${roleRows}</tbody></table></div>

      <h3>The formula</h3>
      <p>For one cell type or sample, let each role contribute the <b>mean normalised expression</b>
      of the genes listed for it (mean = arithmetic mean over those genes; a role with no genes is
      skipped, as described):</p>
      <pre class="mono" style="background:var(--chip);padding:12px;border-radius:6px;overflow-x:auto;font-size:12.5px;line-height:1.7">1  Production
     S_D = sqrt( mean(synthesis) x mean(precursor) )     de novo   (S_D = mean(synthesis) if no precursor genes)
     S_S = sqrt( mean(synthesis) x mean(influx) )        salvage   (only if influx genes exist)
     S   = (S_D + S_S) / 2                               (S = S_D if no salvage term)

2  Net of degradation
     N   = max( 0 , S - mean(degradation) )

3  Availability
     M   = sqrt( N x mean(efflux) )                      (M = N if the ligand has no exporters)

4  Signal strength for one ligand-receptor pair
     Signal = sqrt( M x mean(receptor genes) )           (for a complex: the genes of its subunits)</pre>
      <p>Steps 1, 3 and 4 follow MetaLigand. Step 2 follows the idea in MRCLinkdb of
      producers minus consumers, with negative values set to zero as MetalinksDB does.</p>

      <h3>Worked example: dopamine</h3>
      <p>A hand-picked <b>core gene set</b> is used, because the full dopamine lists on this site are
      noisy (for example 23 degradation genes, including chondroitin-sulfate synthases that come
      from the metabolic model). Core set: synthesis <code>DDC</code>; precursor <code>TH</code>;
      degradation <code>COMT, MAOA, MAOB</code>; efflux <code>SLC18A2</code> (vesicle loading) and
      <code>SLC6A3</code>; receptor <code>DRD2</code>. No salvage genes are used. The expression
      numbers below are <b>invented to show the arithmetic</b>, not real data.</p>
      ${exCells}
      <p>The second cell makes a little dopamine on paper (S = 0.548), but its degradation genes
      are higher than its production, so N is floored at zero and the signal is zero.</p>

      <h3>Try your own numbers</h3>
      <p style="font-size:12.5px;color:var(--muted)">Enter normalised expression values for each
      role, separated by commas (one value per gene). Leave a box empty if the ligand has no genes
      in that role.</p>
      ${box("synthesis", "Synthesis genes")}
      ${box("precursor", "Precursor synthesis genes")}
      ${box("influx", "Influx genes (salvage)")}
      ${box("degradation", "Degradation genes")}
      ${box("efflux", "Efflux genes")}
      ${box("receptor", "Receptor genes")}
      <p><button id="uc-run">Calculate</button>
         <button id="uc-ex1">Example: neuron</button>
         <button id="uc-ex2">Example: other cell</button></p>
      <div id="uc-out"></div>

      <h3>Limits</h3>
      <ul>
        <li>Gene expression is a proxy. Both source papers warn that mRNA and protein levels are
        nonlinearly related, and that metabolite levels also depend on substrates, by-products,
        pH and flux, none of which RNA captures.</li>
        <li>The Efflux role has a built-in assumption: MRCLinkdb's transporter file gives no
        direction, so those rows are stored as efflux and marked with a question mark.</li>
        <li>Source gene lists vary in quality. Transporters and non-enzymes can appear in the
        enzyme lists (MetalinksDB takes them from metabolic models), so a curated core set per
        ligand will usually behave better than the full lists.</li>
        <li>Ions and lipids can cross membranes without transporters, so a missing transporter does
        not prove a ligand is not released.</li>
        <li>This formulation combines two published ones and has not itself been validated.</li>
      </ul>

      <h3>References</h3>
      <ul>
        <li>MetaLigand: prior-knowledge-guided framework for non-peptide ligand mediated cell-cell
        communication. <i>Cell Rep Methods</i> 2025 (PMID 41175863).</li>
        <li>MRCLinkdb: predicting intercellular communication from metabolite-related
        ligand-receptor interactions. <i>BMC Biol</i> 2024 (PMID 38978014).</li>
        <li>MetalinksDB: a resource of metabolite-protein interactions. <i>Brief Bioinform</i>
        2024 (PMID 39038934).</li>
      </ul>`;

    const read = id => document.getElementById("uc-" + id).value
      .split(/[\s,;]+/).filter(Boolean).map(Number).filter(x => isFinite(x) && x >= 0);
    const KEYS = ["synthesis", "precursor", "influx", "degradation", "efflux", "receptor"];
    const run = () => {
      const v = {}; KEYS.forEach(k => v[k] = read(k));
      document.getElementById("uc-out").innerHTML = stepsHtml(score(v));
    };
    const load = cell => { KEYS.forEach(k => document.getElementById("uc-" + k).value = cell[k].join(", ")); run(); };
    const cells = Object.values(EXAMPLE.cells);
    document.getElementById("uc-run").onclick = run;
    document.getElementById("uc-ex1").onclick = () => load(cells[0]);
    document.getElementById("uc-ex2").onclick = () => load(cells[1]);
    run();
  }

  window.buildUsage = build;
  build();
})();
