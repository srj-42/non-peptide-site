/* Gene link-outs for the drawer (connectomeDB-style card): identifiers, Gene Ontology, expression, disease,
   orthologues, aliases. Links only - nothing here enters a score. Data: data/genes.js (build_gene_links.py).
   window.geneCards(symbols, species 'h'|'m'|'any') -> HTML string. */
(function () {
  const esc = s => String(s == null ? "" : s).replace(/[&<>"']/g,
    m => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[m]));
  const enc = encodeURIComponent;
  const a = (href, txt) => `<a href="${esc(href)}" target="_blank" rel="noopener">${esc(txt)}</a>`;
  function card(sym, sp) {
    const G = window.GENES || {};
    const key = sp === "m" ? "m:" + sym : "h:" + sym;
    const g = G[key] || (sp === "any" ? G["m:" + sym] : null);
    if (!g) return "";
    const mouse = !G["h:" + sym] || sp === "m";
    const [up, hgnc, ens, aliases, mgi] = g;
    const L = [];
    if (up) L.push(a("https://www.uniprot.org/uniprotkb/" + up, "UniProt " + up));
    else L.push(a(`https://www.uniprot.org/uniprotkb?query=gene_exact:${enc(sym)}+AND+organism_id:${mouse ? 10090 : 9606}`, "UniProt search"));
    if (hgnc) L.push(a("https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/" + hgnc, hgnc));
    if (ens) L.push(a(`https://www.ensembl.org/${mouse ? "Mus_musculus" : "Homo_sapiens"}/Gene/Summary?g=${ens}`, "Ensembl"));
    const go = up ? a("https://www.ebi.ac.uk/QuickGO/annotations?geneProductId=" + up, "Gene Ontology (QuickGO)")
                  : a("https://www.ebi.ac.uk/QuickGO/search/" + enc(sym), "Gene Ontology (QuickGO)");
    const expr = mouse
      ? (ens ? a("https://www.bgee.org/gene/" + ens, "Bgee") : "")
      : [ens ? a("https://www.proteinatlas.org/" + ens, "Human Protein Atlas") : "", a("https://gtexportal.org/home/gene/" + enc(sym), "GTEx")].filter(Boolean).join(" · ");
    const dis = mouse
      ? a("https://www.alliancegenome.org/search?q=" + enc(sym), "Alliance of Genome Resources")
      : [ens ? a("https://platform.opentargets.org/target/" + ens, "Open Targets") : "", a("https://omim.org/search?search=" + enc(sym), "OMIM")].filter(Boolean).join(" · ");
    const orth = [];
    if (ens) orth.push(a(`https://www.ensembl.org/${mouse ? "Mus_musculus" : "Homo_sapiens"}/Gene/Compara_Ortholog?g=${ens}`, "Ensembl Compara"));
    if (!mouse) orth.push(a("https://www.genenames.org/tools/hcop/#!/?symbol=" + enc(sym), "HGNC HCOP"));
    if (mgi) orth.push(a("https://www.informatics.jax.org/marker/" + mgi, "mouse " + mgi + " (MGI)"));
    if (mouse) orth.push(a("https://www.informatics.jax.org/searchtool/Search.do?query=" + enc(sym), "MGI"));
    const row = (k, v) => v ? `<div><span class="gk">${k}</span> ${v}</div>` : "";
    return `<div class="gcard"><b>${esc(sym)}</b> <span class="gsp">${mouse ? "mouse" : "human"}</span>
      ${row("Identifiers", L.join(" · "))}${row("Gene Ontology", go)}${row("Expression", expr)}${row("Disease", dis)}
      ${row("Orthologues", orth.join(" · "))}${row("Aliases", aliases.length ? esc(aliases.join(", ")) : "")}</div>`;
  }
  window.geneCards = function (syms, sp) {
    const seen = new Set(), html = [];
    for (const s of syms) if (s && !seen.has(s)) { seen.add(s); html.push(card(s, sp)); }
    return html.join("");
  };
  window.geneSection = function (title, syms, sp) {
    const h = window.geneCards(syms, sp);
    return h ? `<dt>${esc(title)}</dt><dd>${h}</dd>` : "";
  };
})();
