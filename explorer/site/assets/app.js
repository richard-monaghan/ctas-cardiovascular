/* ============================================================
   CTAS Explorer — app.js
   ============================================================ */

(function () {
  'use strict';

  // ── State ─────────────────────────────────────────────────

  // ── Query-string fallback ─────────────────────────────────
  // Some static servers use clean-URL redirects (/loci.html?x -> /loci) that drop the query string.
  // Internal links to the loci page therefore also carry the query in the URL fragment, which
  // survives such redirects, and the loci page reads the fragment when the query is empty.
  var QS_PREFIX = '#qs=';
  function withQsFallback(url) {
    var i = url.indexOf('?');
    if (i < 0 || url.indexOf('#') >= 0) return url;
    return url + QS_PREFIX + encodeURIComponent(url.slice(i + 1));
  }
  function currentParams() {
    if (location.search.length > 1) return new URLSearchParams(location.search);
    if (location.hash.indexOf(QS_PREFIX) === 0) {
      try { return new URLSearchParams(decodeURIComponent(location.hash.slice(QS_PREFIX.length))); } catch (e) { /* ignore */ }
    }
    return new URLSearchParams('');
  }
  function clearQsHash(u) { if (u.hash.indexOf(QS_PREFIX) === 0) u.hash = ''; }
  const DATA = { loci: [], domains: [], featured: [], coloc: [], colocPP3: [], colocPP3Map: {}, traitSignals: {}, keywords: {}, memberProfiles: null };

  // Filter state — loci page only
  const FS = {
    query:        '',
    // filterLabel: set when user arrives via ?family= or ?trait= chip nav;
    // provides a semantic display label (e.g. "Coronary artery disease")
    // rather than showing the raw query key in the chips area.
    filterLabel:  null,   // { type: 'family'|'trait', label: string }
    qThreshold:   null,   // null | 'q05' | 'q01' | 'q001' — nested, one at a time
    colocOnly:    false,
    novelOnly:    false,
    locusIds:     null,   // Set of locus_ids when domain/route pre-filter active
    activeFilter: null,   // { type: 'domain'|'route', label: '' }
  };

  let selectedLocusId = null;
  let SCORE_MAX       = 42;   // updated from data

  // Sort state — loci page only.
  // col:  key matching data-sort-col attribute on the active <th>
  // dir:  1 = ascending (A→Z / low→high), -1 = descending (Z→A / high→low)
  const SORT = { col: null, dir: 0 };

  var DEFAULT_LEAD = '91 CTAS loci at genome-wide q<0.05 across 16 Pan-UKBB cardiovascular traits (46 at q<0.01; 21 at q<0.001). Click any row to view locus details.';

  // ── Display-name lookup maps ──────────────────────────────

  var TRAIT_LABELS = {
    'AF_I48':          'Atrial fibrillation',
    'Angina':          'Angina pectoris',
    'Aortic_Aneurysm': 'Aortic aneurysm',
    'Aortic_Valve':    'Aortic valve disease',
    'Atherosclerosis': 'Atherosclerosis',
    'CAD_I25':         'Coronary artery disease',
    'Cardiomyopathy':  'Cardiomyopathy',
    'DVT':             'Phlebitis/thrombophlebitis (ICD-10 I80)',
    'HF_I50':          'Heart failure',
    'Hypertension_401':'Hypertension',
    'Ischemic_Stroke': 'Ischemic stroke',
    'MI':              'Myocardial infarction',
    'Mitral_Valve':    'Mitral valve disease',
    'PAD':             'Other peripheral vascular disease (ICD-10 I73)',
    'TIA':             'Transient ischemic attack',
    'VTE':             'Venous embolism/thrombosis (ICD-10 I82)',
  };

  // Short abbreviations used on the signal-profile bar chart x-axis
  var TRAIT_SHORT = {
    'AF_I48':          'AF',
    'Angina':          'Ang',
    'Aortic_Aneurysm': 'AoAn',
    'Aortic_Valve':    'AoV',
    'Atherosclerosis': 'Ath',
    'CAD_I25':         'CAD',
    'Cardiomyopathy':  'CMP',
    'DVT':             'I80',
    'HF_I50':          'HF',
    'Hypertension_401':'HTN',
    'Ischemic_Stroke': 'IS',
    'MI':              'MI',
    'Mitral_Valve':    'MV',
    'PAD':             'I73',
    'TIA':             'TIA',
    'VTE':             'I82',
  };

  // All 16 trait keys in the fixed canonical Table 1 panel order. The signal
  // profile renders in this order for every locus, so phenotype profiles can be
  // compared position-by-position across loci.
  var TRAIT_KEYS = [
    'CAD_I25','MI','Angina','AF_I48','HF_I50','Cardiomyopathy',
    'Hypertension_401','Atherosclerosis','PAD','Ischemic_Stroke',
    'TIA','DVT','VTE','Aortic_Aneurysm','Aortic_Valve','Mitral_Valve',
  ];

  var TISSUE_LABELS = {
    'aorta':         'Aorta',
    'blood':         'Whole blood',
    'coronary':      'Coronary artery',
    'heart_AA':      'Heart atrial appendage',
    'heart_LV':      'Heart left ventricle',
    'kidney_cortex': 'Kidney cortex',
    'liver':         'Liver',
  };

  // Cardiovascular trait families (eight) — keys match the fam_* source flags
  var FAMILY_LABELS = {
    'coronary':        'Coronary',
    'pressure':        'Pressure',
    'arrhythmia':      'Arrhythmia',
    'heart_failure':   'Heart failure',
    'valve':           'Valve',
    'aortic_vascular': 'Aortic/vascular',
    'cerebrovascular': 'Cerebrovascular',
    'thrombotic':      'Thrombotic',
  };

  // Maps the most-stringent q tier attained to the route IDs used by the
  // clickable q-tier badge in the detail panel.
  var TIER_ROUTE = {
    'q<0.05':  'q05_arch',
    'q<0.01':  'q01_arch',
    'q<0.001': 'q001_core',
  };

  function traitLabel(key) { return TRAIT_LABELS[key]  || key.replace(/_/g, ' '); }
  function tissueLabel(key){ return TISSUE_LABELS[key] || key.replace(/_/g, ' '); }
  function familyLabel(key){ return FAMILY_LABELS[key] || key.replace(/_/g, ' '); }
  // TIER-2 raw primary family, display string ("Broad — no dominant family" when guarded).
  function familyRawDisplay(l){
    if (!l) return '—';
    if (l.family_raw_broad) return 'Broad — no dominant family';
    return l.family_raw ? familyLabel(l.family_raw) : '—';
  }

  // Public, version-neutral locus label: internal id "V2L042" → "Locus 42".
  // The full internal id is retained for routing (data-locus-id) and data joins.
  function publicLocusLabel(id) {
    var m = /(\d+)\s*$/.exec(String(id == null ? '' : id));
    return m ? 'Locus ' + parseInt(m[1], 10) : String(id);
  }

  // FinnGen R12 endpoint code → [short abbreviation, full canonical trait name].
  // Short forms match the trait-chart abbreviations; full names go in title text.
  var FINNGEN_LABELS = {
    'AF':            ['AF',     'Atrial fibrillation'],
    'Angina':        ['Ang',    'Angina pectoris'],
    'CAD':           ['CAD',    'Coronary artery disease'],
    'MI':            ['MI',     'Myocardial infarction'],
    'Hypertension':  ['HTN',    'Hypertension'],
    'HF':            ['HF',     'Heart failure'],
    'VTE':           ['VTE',    'Venous thromboembolism (FinnGen endpoint I9_VTE)'],
    'Stroke':        ['Stroke', 'Stroke'],
    'TIA':           ['TIA',    'Transient ischaemic attack'],
    'Cardiomyopathy':['CMP',    'Cardiomyopathy'],
  };

  // ── Utilities ─────────────────────────────────────────────

  // Cache-busting build tag. Derived from this script's own ?v=<build> query
  // (stamped into the HTML <script src> by scripts/apply_cache_version.py), so a
  // new build changes the HTML, forces a fresh app.js, and in turn versions every
  // JSON fetch below — preventing stale cached data/JS from being served.
  var CTAS_BUILD = (function () {
    try {
      var s = (document.currentScript && document.currentScript.src) || '';
      var m = s.match(/[?&]v=([^&]+)/);
      return m ? m[1] : 'dev';
    } catch (e) { return 'dev'; }
  })();

  async function loadJSON(path) {
    var sep = path.indexOf('?') >= 0 ? '&' : '?';
    const r = await fetch(path + sep + 'v=' + encodeURIComponent(CTAS_BUILD));
    if (!r.ok) throw new Error('Could not load ' + path + ' (' + r.status + ')');
    return r.json();
  }

  function esc(str) {
    return String(str == null ? '' : str)
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;');
  }

  // Maps the q-tier value to the existing badge CSS classes so the established
  // three-colour palette is preserved (most stringent = strongest colour).
  function tierClass(tier) {
    if (tier === 'q<0.001') return 'overlap';
    if (tier === 'q<0.01')  return 'breadth';
    return 'signal';
  }

  function priorLabel(locus) {
    if (locus.not_reported) return { cls: 'novel', text: 'Not in 4 comparators' };
    const n = locus.n_prior_studies || 0;
    if (n >= 3) return { cls: 'known',   text: n + ' prior' };
    if (n === 2) return { cls: 'prior',  text: '2 prior' };
    if (n === 1) return { cls: 'prior',  text: '1 prior' };
    return           { cls: 'known',   text: 'Known' };
  }

  function scorePct(score) {
    return Math.min(100, (score / SCORE_MAX) * 100).toFixed(1);
  }

  // Format a span_kb value for display.
  // span_kb = 0 on single-gene point loci (start === end in the data).
  // We show "< 1 kb" for those rather than "0 kb" or "—".
  function formatSpan(kb) {
    if (kb == null)  return '—';
    if (kb === 0)    return '< 1 kb';
    if (kb >= 1000)  return (kb / 1000).toFixed(2) + ' Mb';
    return kb.toFixed(1) + ' kb';
  }

  // PP4 display: 3 dp normally; 4 dp when 3 dp would round to 1.000 (e.g. 0.9999)
  function pp4Display(v) {
    if (v == null) return '—';
    return v >= 0.9995 ? v.toFixed(4) : v.toFixed(3);
  }

  // Span chip (non-clickable)
  function chip(text, cls) {
    return '<span class="chip ' + (cls || '') + '">' + esc(text) + '</span>';
  }

  // Anchor chip (clickable)
  function chipLink(text, cls, href) {
    return '<a class="chip chip-clickable ' + (cls || '') + '" href="' + esc(withQsFallback(href)) + '">' + esc(text) + '</a>';
  }

  // Span badge (non-clickable)
  function badge(text, cls) {
    return '<span class="badge ' + (cls || '') + '">' + esc(text) + '</span>';
  }

  // Anchor badge (clickable)
  function badgeLink(text, cls, href) {
    return '<a class="badge badge-clickable ' + (cls || '') + '" href="' + esc(withQsFallback(href)) + '">' + esc(text) + '</a>';
  }

  // ── Contextual help: one reusable, accessible popover ─────
  // Central help-text dictionary (single source of truth; the Method page
  // remains the scientific authority — these strings match its definitions).
  var HELP_TEXT = {
    ctas: 'The CTAS statistic combines residual signal magnitude with effective breadth: CTAS = M × log2(1 + PR). Higher values reflect stronger and/or more widely distributed cross-trait evidence after dependence correction. BH q reports genome-wide significance.',
    magnitude: 'The total amount of positive residual evidence after decorrelating the complete 16-trait profile. Larger values indicate stronger aggregate residual signal.',
    pr: 'The effective number of decorrelated, or whitened, dimensions carrying positive evidence. Larger values indicate that the cross-trait signal is distributed more broadly.',
    family_neff: 'The effective breadth of the descriptive profile across the eight named cardiovascular families. Larger values indicate a more evenly distributed family profile.',
    leading_family: 'The family or co-leading families contributing the largest shares of the size-adjusted descriptive family profile. Two families are shown when the second contributes at least 90% as much as the first. This summarises the full profile and is not a causal classification.',
    stability: 'The number of 16 leave-one-trait-out sensitivity runs in which the locus remained prioritised at q<0.05. A value of 16/16 means the locus was retained whichever single trait was removed. These precomputed runs use the raw-Spearman sensitivity matrix and a 10⁶-draw sensitivity null; they were not rerun under the primary analysis.',
    traits_p01: 'Number of the 16 original traits with raw MAGMA gene-level P<0.01 for the representative gene. This provides a compact summary of the named-trait profile.',
    traits_p001: 'Number of the 16 original traits with raw MAGMA gene-level P<0.001. This provides the same named-trait summary at a stricter threshold.',
    families_p01: 'Number of the eight cardiovascular families containing at least one trait with raw gene-level P<0.01.',
    span: 'Genomic width of the locked locus interval from its start to end coordinate.',
    signal_profile: 'The representative gene’s raw MAGMA −log10P profile across all 16 traits in a fixed order. CTAS transforms and dependence-corrects this complete profile jointly.',
    representative: 'The highest-CTAS significant gene used as the reproducible label for the locus. It is the statistical representative used throughout the Explorer and manuscript.',
    qtier: 'The most stringent genome-wide Benjamini–Hochberg q threshold reached by the locus. The three nested levels provide progressively stricter views of the prioritised architecture.',
    prior: 'Number of the four prespecified comparator studies—Qiao, Liu, Foley (coordinate containment) and PLEIO—in which the locus was matched, each on its native unit. This summarises convergence with earlier cross-trait cardiovascular studies; absence is not evidence of biological novelty.',
    coloc: 'The strongest displayed colocalisation result. PP4 supports a shared-variant model, while PP3 indicates association and expression signals arising from different local variants.',
    landscape: 'Each point represents one prioritised locus, positioned in genomic order and plotted by its representative gene’s CTAS. Colour shows its most stringent BH q tier; selecting a point opens the locus browser.',
    family_browse: 'A locus appears in a family when its representative gene has at least one original trait in that family with raw gene-level P<0.01.',
    coloc_section: 'PP4 measures support for the cardiovascular and eQTL associations sharing a causal variant. PP3 measures support for both signals being present through different local variants.',
    finngen: 'Representative-gene FinnGen support: the number of six clinical families (coronary, arrhythmia, myocardial, pressure, cerebrovascular, venous) in which at least one FinnGen endpoint has gene-level P<0.01 for the statistical locus representative, shown together with the number of the ten endpoints. Endpoints are grouped so that correlated endpoints are not counted as independent support. This summarises external support in an independent cohort. “Not evaluable” means that no usable gene-level result was available; it is not scored as zero.',
    rep_evidence: 'CTAS, magnitude and effective breadth for the highest-CTAS gene used as the statistical locus representative. Multi-gene loci also expose the individual CTAS components and raw 16-trait profiles of every significant member gene below; these gene-level statistics are not combined into a locus score.',
    family_raw: 'The cardiovascular family carrying the strongest evidence in the representative gene’s original MAGMA profile, before dependence correction. This is a descriptive summary of the raw trait profile, not a causal classification or the CTAS selection rule. Loci resting on a single nominally significant trait, or with no family reaching about 30% of the raw evidence, are shown as “Broad — no dominant family” rather than forced into a family.',
    known_family: 'The cardiovascular context previously associated with this gene or genomic region in the curated literature. It is contextual evidence rather than a CTAS-derived result, and it does not assign a causal gene or mechanism.',
    family_ctas_weighted: 'The family carrying the largest share of the representative gene’s positive residual evidence after the fitted cross-trait dependence has been discounted. It can differ from the raw-profile family because correlated signals are not counted repeatedly. It is a descriptive projection of the whitened profile, not a separate significance test.'
  };
  // "Top family" (table) now reports the raw primary family.
  HELP_TEXT.top_family = HELP_TEXT.family_raw;
  var HELP_LABEL = {
    ctas: 'CTAS', magnitude: 'magnitude', pr: 'effective breadth', family_neff: 'family Neff',
    leading_family: 'leading family profile', top_family: 'the raw-profile leading family', stability: 'stability',
    family_raw: 'the raw-profile leading family', known_family: 'prior cardiovascular context', family_ctas_weighted: 'the leading family after dependence correction',
    traits_p01: 'traits at P<0.01', traits_p001: 'traits at P<0.001', families_p01: 'families at P<0.01',
    span: 'span', signal_profile: 'the signal profile', representative: 'representative gene',
    qtier: 'q tier', prior: 'prior support', coloc: 'colocalisation', landscape: 'the locus landscape',
    family_browse: 'browse by disease family', coloc_section: 'colocalisation',
    finngen: 'representative-gene FinnGen support', rep_evidence: 'representative-gene CTAS evidence'
  };

  // Returns the markup for one help trigger. `extraCls` e.g. 'help-inline'.
  function helpBtn(key, extraCls) {
    var label = HELP_LABEL[key] || key;
    return '<button type="button" class="help-btn' + (extraCls ? ' ' + extraCls : '') + '"' +
           ' data-help="' + key + '" aria-label="Explain ' + esc(label) + '"' +
           ' aria-expanded="false">?</button>';
  }

  var _helpPop = null, _helpBtnOpen = null, _helpPinned = false;
  var _helpHoverTimer = null;            // pending 400 ms hover-open timer
  var HELP_HOVER_DELAY = 400;            // ms; pointer hover only (click/focus are immediate)
  function _cancelHoverOpen() { if (_helpHoverTimer) { clearTimeout(_helpHoverTimer); _helpHoverTimer = null; } }

  function _helpEl() {
    if (_helpPop) return _helpPop;
    _helpPop = document.createElement('div');
    _helpPop.className = 'help-popover';
    _helpPop.id = 'help-popover';
    _helpPop.setAttribute('role', 'tooltip');
    document.body.appendChild(_helpPop);
    return _helpPop;
  }
  function _positionHelp(btn) {
    var pop = _helpEl();
    var r = btn.getBoundingClientRect();
    var pw = pop.offsetWidth, ph = pop.offsetHeight, m = 8;
    var left = r.left + r.width / 2 - pw / 2;
    left = Math.max(m, Math.min(left, window.innerWidth - pw - m));
    var top = r.bottom + 6;
    if (top + ph > window.innerHeight - m) {
      var above = r.top - ph - 6;
      top = above >= m ? above : Math.max(m, window.innerHeight - ph - m);
    }
    pop.style.left = Math.round(left) + 'px';
    pop.style.top = Math.round(top) + 'px';
  }
  function openHelp(btn, pinned) {
    var key = btn.getAttribute('data-help');
    var txt = HELP_TEXT[key];
    if (!txt) return;
    if (_helpBtnOpen && _helpBtnOpen !== btn) {
      _helpBtnOpen.setAttribute('aria-expanded', 'false');
      _helpBtnOpen.removeAttribute('aria-describedby');
    }
    var pop = _helpEl();
    pop.textContent = txt;
    _helpBtnOpen = btn;
    _helpPinned = !!pinned;
    btn.setAttribute('aria-expanded', 'true');
    btn.setAttribute('aria-describedby', 'help-popover');
    pop.classList.add('open');
    _positionHelp(btn);
  }
  function closeHelp() {
    _cancelHoverOpen();
    if (_helpPop) _helpPop.classList.remove('open');
    if (_helpBtnOpen) {
      _helpBtnOpen.setAttribute('aria-expanded', 'false');
      _helpBtnOpen.removeAttribute('aria-describedby');
    }
    _helpBtnOpen = null;
    _helpPinned = false;
  }
  function _closestHelp(t) { return t && t.closest ? t.closest('.help-btn') : null; }

  // Click / tap: capture phase + stopPropagation so a help button inside a
  // sortable header or a locus row never triggers sorting or the drawer.
  document.addEventListener('click', function (e) {
    var btn = _closestHelp(e.target);
    if (btn) {
      e.preventDefault();
      e.stopPropagation();
      _cancelHoverOpen();             // click is immediate; drop any pending hover-open
      var pointer = e.detail !== 0;   // detail===0 => keyboard-activated
      if (pointer) {
        if (_helpBtnOpen === btn && _helpPinned) closeHelp();
        else openHelp(btn, true);
      } else {
        openHelp(btn, false);         // keyboard: transient, tied to focus
      }
      return;
    }
    if (_helpBtnOpen) closeHelp();     // click-away
  }, true);

  document.addEventListener('keydown', function (e) {
    if ((e.key === 'Enter' || e.key === ' ' || e.key === 'Spacebar') && _closestHelp(e.target)) {
      e.stopPropagation();             // block header sort on Enter/Space; native click still fires
    }
    if (e.key === 'Escape' && _helpBtnOpen) {
      var b = _helpBtnOpen;
      closeHelp();
      if (b && b.focus) b.focus();
    }
  }, true);

  // Pointer hover: open only after HELP_HOVER_DELAY of continuous hover on ONE control.
  // Any new hover / leave / click / focus / scroll cancels a still-pending timer, so
  // moving across controls never opens the wrong popover and never flickers.
  document.addEventListener('mouseover', function (e) {
    var btn = _closestHelp(e.target);
    if (!btn || _helpPinned) return;
    _cancelHoverOpen();                        // reset delay for the newly hovered control
    _helpHoverTimer = setTimeout(function () {
      _helpHoverTimer = null;
      if (!_helpPinned && document.contains(btn)) openHelp(btn, false);
    }, HELP_HOVER_DELAY);
  });
  document.addEventListener('mouseout', function (e) {
    var btn = _closestHelp(e.target);
    if (!btn) return;
    _cancelHoverOpen();                        // left before the delay elapsed → no open
    if (!_helpPinned && _helpBtnOpen === btn) closeHelp();
  });
  document.addEventListener('focusin', function (e) {
    var btn = _closestHelp(e.target);
    if (btn) { _cancelHoverOpen(); openHelp(btn, false); }   // keyboard focus is immediate
  });
  document.addEventListener('focusout', function (e) {
    var btn = _closestHelp(e.target);
    if (btn && !_helpPinned && _helpBtnOpen === btn) closeHelp();
  });
  window.addEventListener('scroll', function () {
    _cancelHoverOpen();                        // scrolling away must not open a stale popover
    if (!_helpBtnOpen) return;
    if (_helpPinned) _positionHelp(_helpBtnOpen);
    else closeHelp();
  }, true);
  window.addEventListener('resize', function () { if (_helpBtnOpen) _positionHelp(_helpBtnOpen); });

  // ── Nav title / page header helpers ──────────────────────

  function setNavTitle(title) {
    var el = document.getElementById('nav-title');
    if (el) el.textContent = title;
    document.title = (title === 'CTAS Explorer')
      ? 'CTAS Explorer'
      : title + ' — CTAS Explorer';
  }

  function setPageHeader(label, title, lead) {
    var elLabel = document.getElementById('page-label');
    var elTitle = document.getElementById('page-title');
    var elLead  = document.getElementById('page-lead');
    if (elLabel && label !== undefined) elLabel.textContent = label;
    if (elTitle && title !== undefined) elTitle.textContent = title;
    if (elLead  && lead  !== undefined) elLead.textContent  = lead;
  }

  function resetPageHeader() {
    setNavTitle('Loci');
    setPageHeader('Loci browser', 'Cross-trait loci', DEFAULT_LEAD);
  }

  // ── Data normalization ────────────────────────────────────
  function normalizeLoci(loci) {
    return loci.map(function (l) {
      l.genes            = [].concat(l.genes            || []);
      l.family_list_p01  = [].concat(l.family_list_p01  || []);
      l.family_list_p001 = [].concat(l.family_list_p001 || []);
      l.trait_list_p01   = [].concat(l.trait_list_p01   || []);
      l.trait_list_p001  = [].concat(l.trait_list_p001  || []);
      l.coloc_genes      = [].concat(l.coloc_genes      || []);
      return l;
    });
  }

  // Join keyword annotations onto locus objects by locus_id.
  function joinKeywords(loci, keywords) {
    loci.forEach(function (l) {
      l._kw = keywords[String(l.locus_id)] || null;
    });
  }

  // Return the keyword annotation for a locus (or null).
  function kwFor(locus) { return locus._kw || null; }

  // Build intuitive free-text search terms for the coloc layer.
  // This keeps the browser searchable by user-facing descriptors
  // (e.g. "PP4", "regional", tissues, traits) rather than raw locus IDs.
  function colocSearchTerms(locus) {
    var terms = [];
    var pp3Entry = DATA.colocPP3Map[locus.locus_id];

    if (locus.has_coloc) {
      terms.push('coloc', 'colocalization', 'colocalisation', 'pp4', 'shared signal');
      if (locus.locus_best_tissue) terms.push(tissueLabel(locus.locus_best_tissue), locus.locus_best_tissue);
      if (locus.locus_best_trait)  terms.push(traitLabel(locus.locus_best_trait), locus.locus_best_trait);
      (locus.coloc_genes || []).forEach(function (g) { terms.push(g); });
    }

    if (pp3Entry && pp3Entry.best_pp3 != null && pp3Entry.best_pp3 >= 0.5) {
      terms.push('pp3', 'regional', 'regional overlap', 'distinct local variants');
      if (pp3Entry.best_pp3_tissue) terms.push(tissueLabel(pp3Entry.best_pp3_tissue), pp3Entry.best_pp3_tissue);
      if (pp3Entry.best_pp3_trait)  terms.push(traitLabel(pp3Entry.best_pp3_trait), pp3Entry.best_pp3_trait);
      if (pp3Entry.best_pp3_gene)   terms.push(pp3Entry.best_pp3_gene);
    }

    if (locus.coloc_class === 'not_assessed') {
      terms.push('coloc not assessed', 'not assessed');
    } else if (!locus.has_coloc && (!pp3Entry || pp3Entry.coloc_class !== 'PP3-only')) {
      terms.push('no coloc', 'no colocalization', 'no colocalisation');
    }

    return terms;
  }

  // Free-text terms for the Stage-2 evidence layer (FinnGen, comparators,
  // biological-context gene) so the drawer's evidence is reachable by search.
  function evidenceSearchTerms(locus) {
    var t = [];
    if (locus.biological_context_gene) t.push(locus.biological_context_gene);
    if (locus.finngen_sig_traits) t.push(String(locus.finngen_sig_traits).replace(/;/g, ' '));
    if (locus.finngen_available) t.push('finngen');
    if (locus.in_qiao_2025)  t.push('Qiao');
    if (locus.in_liu_2023)   t.push('Liu');
    if (locus.in_foley_2021) t.push('Foley');
    if (locus.in_pleio)      t.push('PLEIO');
    return t;
  }

  // Classify a display_keyword as gene-scope ('gene') or locus-context ('locus')
  // by checking whether any word token (length ≥ 3) of the display keyword also
  // appears in the gene_keywords array.
  function classifyDisplayKeyword(dkw, geneKws) {
    var dTokens = dkw.toLowerCase().split(/[^a-z0-9]+/).filter(function (w) { return w.length >= 3; });
    var gTokenSet = new Set();
    (geneKws || []).forEach(function (gkw) {
      gkw.toLowerCase().split(/[^a-z0-9]+/).forEach(function (w) {
        if (w.length >= 3) gTokenSet.add(w);
      });
    });
    return dTokens.some(function (t) { return gTokenSet.has(t); }) ? 'gene' : 'locus';
  }

  // Render 3 keyword chips for the detail panel header row.
  function renderKwHeaderChips(kw) {
    if (!kw || !(kw.display_keywords || []).length) return '';
    var geneKws = kw.gene_keywords || [];
    var chips = kw.display_keywords.map(function (dkw) {
      var scope = classifyDisplayKeyword(dkw, geneKws);
      return chip(dkw, scope === 'gene' ? 'chip-kw-gene' : 'chip-kw-locus');
    });
    return '<div class="detail-kw-header-row">' + chips.join('') + '</div>';
  }

  // Render the full keyword section for the detail panel (Lead gene / Locus context / Biological context).
  function renderKeywordSection(kw) {
    if (!kw) return '';
    var geneKws  = kw.gene_keywords || [];
    var locusKws = kw.locus_context_keywords || [];
    var bioCtxs  = kw.biological_context || [];

    var geneGroup = geneKws.length
      ? '<div class="kw-group"><span class="kw-group-label">Lead gene</span>' +
        geneKws.map(function (k) { return chip(k, 'chip-kw-gene'); }).join('') + '</div>'
      : '';

    var locusGroup = locusKws.length
      ? '<div class="kw-group"><span class="kw-group-label">Locus context</span>' +
        locusKws.map(function (k) { return chip(k, 'chip-kw-locus'); }).join('') + '</div>'
      : '';

    var bioGroup = bioCtxs.length
      ? '<div class="kw-group"><span class="kw-group-label">Biological context</span>' +
        bioCtxs.map(function (k) { return chip(k, 'chip-bio-ctx'); }).join('') + '</div>'
      : '';

    var candidateHtml = '';
    if (kw.locus_candidate && kw.locus_context_note) {
      candidateHtml = [
        '<div class="kw-candidate-box">',
        '  <span class="kw-candidate-gene-label">Candidate gene: </span>',
        '  <span class="kw-candidate-gene">' + esc(kw.locus_candidate) + '</span>',
        '  <p class="kw-candidate-note-text">' + esc(kw.locus_context_note) + '</p>',
        '</div>',
      ].join('\n');
    }

    if (!geneGroup && !locusGroup && !bioGroup) return '';

    return [
      '<div class="detail-section">',
      '  <div class="detail-section-title">Keywords</div>',
      geneGroup,
      locusGroup,
      bioGroup,
      candidateHtml,
      '</div>',
    ].join('\n');
  }

  // ── Navigation ────────────────────────────────────────────

  function setActiveNav() {
    var slug   = location.pathname.split('/').pop().replace(/\.html$/, '') || 'index';
    var isHome = (slug === '' || slug === 'index');

    var backBtn = document.getElementById('nav-back');
    if (backBtn) backBtn.hidden = isHome;

    if (isHome)           setNavTitle('CTAS Explorer');
    else if (slug==='loci')   setNavTitle('Loci');
    else if (slug==='method') setNavTitle('Method');

    document.querySelectorAll('.nav-link').forEach(function (a) {
      var href = a.getAttribute('href').replace(/\.html$/, '');
      a.classList.toggle('active', href === slug);
    });
  }

  function initBackButton() {
    var backBtn = document.getElementById('nav-back');
    if (!backBtn) return;
    backBtn.addEventListener('click', function () {
      if (history.length > 1 && document.referrer &&
          document.referrer.indexOf(location.hostname) !== -1) {
        history.back();
      } else {
        location.href = 'index.html';
      }
    });
  }

  function initMobileNav() {
    var toggle = document.querySelector('.nav-toggle');
    var links  = document.querySelector('.nav-links');
    if (!toggle || !links) return;
    toggle.addEventListener('click', function () {
      var open = links.classList.toggle('open');
      toggle.setAttribute('aria-expanded', String(open));
    });
    document.addEventListener('click', function (e) {
      if (!toggle.contains(e.target) && !links.contains(e.target)) {
        links.classList.remove('open');
        toggle.setAttribute('aria-expanded', 'false');
      }
    });
  }

  // ── Locus Landscape Plot (index.html) ────────────────────
  // Full-width SVG lollipop chart: x = genomic order, y = CTAS,
  // colour = tier.  Hover shows a tooltip; click navigates to the loci
  // browser with that locus's detail panel auto-opened via ?locus=<id>.

  function renderLandscapePlot(loci) {
    var wrap = document.getElementById('landscape-wrap');
    if (!wrap || !loci.length) return;

    // ── Genomic sort ──────────────────────────────────────
    function chrNum(c) {
      var n = parseInt(c, 10);
      return isNaN(n) ? (c === 'X' ? 23 : c === 'Y' ? 24 : 25) : n;
    }
    var sorted = loci.slice().sort(function (a, b) {
      var ca = chrNum(a.chr), cb = chrNum(b.chr);
      if (ca !== cb) return ca - cb;
      return (a.start || 0) - (b.start || 0);
    });

    // ── Layout ────────────────────────────────────────────
    // Fixed viewBox; CSS scales SVG to container width.
    var VW = 1060, VH = 210;
    var mL = 46, mR = 14, mT = 18, mB = 44;
    var cW = VW - mL - mR;     // chart width  (1000)
    var cH = VH - mT - mB;     // chart height (148)
    var n  = sorted.length;
    var maxScore = sorted.reduce(function (m, l) { return Math.max(m, l.ctas_score); }, 0);
    var yMax = maxScore * 1.06;  // 6 % headroom at the top
    var ptW  = cW / n;           // x-slot width per locus

    function xPos(i) { return mL + (i + 0.5) * ptW; }
    function yPos(s) { return mT + cH - (s / yMax) * cH; }

    // ── Point array ───────────────────────────────────────
    var pts = sorted.map(function (locus, i) {
      return { locus: locus, x: xPos(i), y: yPos(locus.ctas_score) };
    });

    // ── Chromosome groups ─────────────────────────────────
    var chrGroups = [];
    var curChr = null, chrStart = 0;
    pts.forEach(function (pt, i) {
      if (pt.locus.chr !== curChr) {
        if (curChr !== null) chrGroups.push({ chr: curChr, from: chrStart, to: i - 1 });
        curChr = pt.locus.chr;
        chrStart = i;
      }
    });
    chrGroups.push({ chr: curChr, from: chrStart, to: n - 1 });

    // ── q-tier palette (most stringent = strongest colour) ─
    var FILL = {
      'q<0.001': '#4A7A58',
      'q<0.01':  '#4A7FA5',
      'q<0.05':  '#7E8A9A',
    };

    // ── Y-axis ticks ──────────────────────────────────────
    var step = maxScore > 30 ? 10 : 5;
    var yTicks = [];
    for (var v = 0; v <= Math.ceil(maxScore / step) * step; v += step) yTicks.push(v);

    // ── Build SVG string ─────────────────────────────────
    var s = [];
    s.push('<svg class="landscape-svg" viewBox="0 0 ' + VW + ' ' + VH + '"' +
           ' xmlns="http://www.w3.org/2000/svg">');

    // Alternating chromosome band shading (even-indexed groups)
    chrGroups.forEach(function (g, gi) {
      if (gi % 2 !== 1) return;
      var x1 = pts[g.from].x - ptW * 0.5;
      var x2 = pts[g.to].x   + ptW * 0.5;
      s.push('<rect x="' + x1.toFixed(1) + '" y="' + mT + '"' +
             ' width="' + (x2 - x1).toFixed(1) + '" height="' + cH + '"' +
             ' fill="#EDEBE4" opacity="0.55"/>');
    });

    // Chromosome boundary lines + labels
    chrGroups.forEach(function (g, gi) {
      if (gi > 0) {
        var bx = (pts[g.from].x - ptW * 0.5).toFixed(1);
        s.push('<line x1="' + bx + '" y1="' + mT + '" x2="' + bx + '" y2="' + (mT + cH) +
               '" stroke="var(--border)" stroke-width="0.75" opacity="0.5"/>');
      }
      // Label: only when the group is wide enough (≥ 2 loci)
      if (g.to > g.from) {
        var cx = ((pts[g.from].x + pts[g.to].x) / 2).toFixed(1);
        s.push('<text class="ls-chrlab" x="' + cx + '" y="' + (mT + cH + 18) +
               '" text-anchor="middle">' + esc(g.chr) + '</text>');
      }
    });

    // X-axis descriptor (bottom)
    s.push('<text class="ls-xlab" x="' + (mL + cW / 2).toFixed(1) + '"' +
           ' y="' + (VH - 4) + '" text-anchor="middle">Genomic order (chromosomal)</text>');

    // Y gridlines and labels
    yTicks.forEach(function (v) {
      var y = yPos(v).toFixed(1);
      s.push('<line class="ls-grid" x1="' + mL + '" y1="' + y + '" x2="' + (mL + cW) + '" y2="' + y + '"/>');
      s.push('<text class="ls-ylab" x="' + (mL - 5) + '" y="' + (yPos(v) + 3.5).toFixed(1) + '"' +
             ' text-anchor="end">' + v + '</text>');
    });

    // Axes
    s.push('<line class="ls-axis" x1="' + mL + '" y1="' + (mT + cH) + '"' +
           ' x2="' + (mL + cW) + '" y2="' + (mT + cH) + '"/>');
    s.push('<line class="ls-axis" x1="' + mL + '" y1="' + mT + '"' +
           ' x2="' + mL + '" y2="' + (mT + cH) + '"/>');

    // Y-axis title (rotated)
    var ytX = mL - 34, ytY = (mT + cH / 2).toFixed(1);
    s.push('<text class="ls-ytitle"' +
           ' transform="rotate(-90,' + ytX + ',' + ytY + ')"' +
           ' x="' + ytX + '" y="' + ytY + '" text-anchor="middle">CTAS</text>');

    // Lollipop stems (rendered before dots)
    pts.forEach(function (pt) {
      var fill = FILL[pt.locus.tier] || '#7E8A9A';
      s.push('<line class="ls-stem"' +
             ' x1="' + pt.x.toFixed(1) + '" y1="' + pt.y.toFixed(1) + '"' +
             ' x2="' + pt.x.toFixed(1) + '" y2="' + (mT + cH) + '"' +
             ' stroke="' + fill + '"/>');
    });

    // Dots — q-tier z-order: least stringent behind most stringent
    ['q<0.05', 'q<0.01', 'q<0.001'].forEach(function (tier) {
      pts.filter(function (pt) { return pt.locus.tier === tier; })
         .forEach(function (pt) {
           var fill = FILL[tier] || '#7E8A9A';
           s.push('<circle class="ls-dot"' +
                  ' data-lid="' + pt.locus.locus_id + '"' +
                  ' cx="' + pt.x.toFixed(1) + '" cy="' + pt.y.toFixed(1) + '"' +
                  ' r="4.5" fill="' + fill + '"/>');
         });
    });

    s.push('</svg>');
    wrap.innerHTML = s.join('');

    // ── Legend ────────────────────────────────────────────
    var legend = document.getElementById('landscape-legend');
    if (legend) {
      legend.innerHTML = ['q<0.001', 'q<0.01', 'q<0.05'].map(function (tier) {
        return '<span class="ls-leg">' +
          '<svg width="10" height="10" viewBox="0 0 10 10">' +
          '<circle cx="5" cy="5" r="4.5" fill="' + (FILL[tier] || '#999') + '"/>' +
          '</svg>' + esc(tier) + '</span>';
      }).join('');
    }

    // ── Tooltip element ───────────────────────────────────
    var tipEl = document.createElement('div');
    tipEl.className = 'ls-tooltip';
    tipEl.hidden = true;
    document.body.appendChild(tipEl);

    // ── Interaction ───────────────────────────────────────
    var activeDot = null;

    function clearHover() {
      if (activeDot) { activeDot.setAttribute('r', '4.5'); activeDot = null; }
      tipEl.hidden = true;
    }

    wrap.addEventListener('mousemove', function (e) {
      var dot = e.target.closest ? e.target.closest('.ls-dot') : null;
      if (!dot) { clearHover(); return; }

      var id    = dot.getAttribute('data-lid');
      var locus = DATA.loci.find(function (l) { return l.locus_id === id; });
      if (!locus) { clearHover(); return; }

      // Enlarge dot on hover
      if (activeDot !== dot) {
        if (activeDot) activeDot.setAttribute('r', '4.5');
        dot.setAttribute('r', '6.5');
        activeDot = dot;
      }

      // Populate tooltip
      var fill = FILL[locus.tier] || '#7E8A9A';
      tipEl.innerHTML =
        '<div class="ls-tip-gene">' + esc(locus.lead_gene) + '</div>' +
        '<div class="ls-tip-meta">' +
          '<span class="ls-tip-swatch" style="background:' + fill + '"></span>' +
          '<span>' + esc(locus.tier) + '</span>' +
        '</div>' +
        '<div class="ls-tip-score">CTAS&thinsp;<strong>' + locus.ctas_score.toFixed(3) + '</strong></div>';

      // Position near cursor, keeping within viewport
      var tx = e.clientX + 16;
      var ty = e.clientY - 14;
      if (tx + 200 > window.innerWidth)  tx = e.clientX - 200 - 10;
      if (ty + 90  > window.innerHeight) ty = e.clientY - 90;
      tipEl.style.left = tx + 'px';
      tipEl.style.top  = ty + 'px';
      tipEl.hidden = false;
    });

    wrap.addEventListener('mouseleave', clearHover);

    // Click → loci page, auto-open that locus detail.
    // Store the target id in sessionStorage before navigating so the loci page
    // can read it synchronously after the table renders — no URL param parsing,
    // no scroll-restoration race, no rAF timing issues.
    wrap.addEventListener('click', function (e) {
      var dot = e.target.closest ? e.target.closest('.ls-dot') : null;
      if (!dot) return;
      sessionStorage.setItem('ctas_open_locus', dot.getAttribute('data-lid'));
      location.href = 'loci.html';
    });
  }

  // ── Domain Cards (index.html) ─────────────────────────────

  function renderDomainCards(domains) {
    var grid = document.getElementById('domain-grid');
    if (!grid) return;

    var usedGenes = new Set();

    grid.innerHTML = domains.map(function (d) {
      var genes = [];
      var lids  = [].concat(d.locus_ids || []);
      for (var i = 0; i < lids.length && genes.length < 3; i++) {
        var locus = DATA.loci.find(function (l) { return l.locus_id === lids[i]; });
        if (locus && !usedGenes.has(locus.lead_gene)) {
          genes.push(locus.lead_gene);
          usedGenes.add(locus.lead_gene);
        }
      }
      var examples = genes.join(', ') || '—';
      return [
        '<article class="domain-card" role="listitem" tabindex="0"',
        '         data-domain-id="' + esc(d.id) + '">',
        '  <div class="domain-card-name">' + esc(d.label) + '</div>',
        '  <div class="domain-card-count">' + d.n_loci + ' loci</div>',
        '  <div class="domain-card-genes">' + esc(examples) + '</div>',
        '</article>',
      ].join('\n');
    }).join('');

    grid.addEventListener('click', function (e) {
      var card = e.target.closest('[data-domain-id]');
      if (card) location.href = withQsFallback('loci.html?domain=' + encodeURIComponent(card.dataset.domainId));
    });
    grid.addEventListener('keydown', function (e) {
      if (e.key === 'Enter' || e.key === ' ') {
        var card = e.target.closest('[data-domain-id]');
        if (card) location.href = withQsFallback('loci.html?domain=' + encodeURIComponent(card.dataset.domainId));
      }
    });
  }

  // ── Featured Route Cards (index.html) ────────────────────

  function renderFeaturedCards(featured) {
    var grid = document.getElementById('routes-grid');
    if (!grid) return;

    grid.innerHTML = featured.map(function (r) {
      return [
        '<article class="route-card" role="listitem" tabindex="0"',
        '         data-route-id="' + esc(r.id) + '">',
        '  <span class="route-tag">' + esc(r.tag) + '</span>',
        '  <div class="route-title">' + esc(r.label) + '</div>',
        '  <p class="route-desc">' + esc(r.description) + '</p>',
        '  <div class="route-footer">',
        '    <span class="route-count">' + r.n_loci + ' loci</span>',
        '    <span class="route-link" aria-hidden="true">View &#8594;</span>',
        '  </div>',
        '</article>',
      ].join('\n');
    }).join('');

    grid.addEventListener('click', function (e) {
      var card = e.target.closest('[data-route-id]');
      if (card) location.href = withQsFallback('loci.html?route=' + encodeURIComponent(card.dataset.routeId));
    });
    grid.addEventListener('keydown', function (e) {
      if (e.key === 'Enter' || e.key === ' ') {
        var card = e.target.closest('[data-route-id]');
        if (card) location.href = withQsFallback('loci.html?route=' + encodeURIComponent(card.dataset.routeId));
      }
    });
  }

  // ── Hero Search (index.html) ──────────────────────────────

  function initHeroSearch(loci) {
    var input   = document.querySelector('.search-input');
    var note    = document.querySelector('.search-note');
    var results = document.getElementById('search-results');
    if (!input || !results) return;

    if (note) note.textContent = 'Searches genes, traits, biology, and coloc descriptors.';

    function doSearch(q) {
      if (!q || q.length < 2) { results.hidden = true; return; }
      var ql = q.toLowerCase();
      var hits = loci.filter(function (l) {
        var kw = kwFor(l);
        var colocTerms = colocSearchTerms(l);
        return l.lead_gene.toLowerCase().includes(ql)
          || l.genes.some(function (g) { return g.toLowerCase().includes(ql); })
          || (l.annotation_note && l.annotation_note.toLowerCase().includes(ql))
          || (l.family_list_p01 || []).some(function (f) {
               return f.toLowerCase().includes(ql) || familyLabel(f).toLowerCase().includes(ql);
             })
          || (l.leading_family_keys || []).some(function (k) {
               return k.toLowerCase().includes(ql) || familyLabel(k).toLowerCase().includes(ql);
             })
          || (l.trait_list_p01 || []).some(function (t) {
               return t.toLowerCase().includes(ql) || traitLabel(t).toLowerCase().includes(ql);
             })
          || (kw && (kw.search_keywords || []).some(function (k) { return k.toLowerCase().includes(ql); }))
          || (kw && (kw.biological_context || []).some(function (k) { return k.toLowerCase().includes(ql); }))
          || (kw && kw.locus_candidate && kw.locus_candidate.toLowerCase().includes(ql))
          || colocTerms.some(function (t) { return String(t).toLowerCase().includes(ql); })
          || evidenceSearchTerms(l).some(function (t) { return t.toLowerCase().includes(ql); });
      }).slice(0, 8);

      if (!hits.length) { results.hidden = true; return; }
      results.hidden = false;
      results.innerHTML = hits.map(function (l) {
        var tc  = tierClass(l.tier);
        var kw  = kwFor(l);
        var snippet = (kw && (kw.display_keywords || []).length)
          ? kw.display_keywords.join(' · ')
          : (l.annotation_note || '').slice(0, 64) + ((l.annotation_note || '').length > 64 ? '…' : '');
        return [
          '<a class="search-result"',
          '   href="' + esc(withQsFallback('loci.html?q=' + encodeURIComponent(l.lead_gene))) + '">',
          '  <span class="search-result-gene">' + esc(l.lead_gene) + '</span>',
          '  <span class="search-result-note">' + esc(snippet) + '</span>',
          '  <span class="tier ' + tc + '">' + esc(l.tier) + '</span>',
          '</a>',
        ].join('\n');
      }).join('');
    }

    input.addEventListener('input', function () { doSearch(input.value.trim()); });
    input.addEventListener('keydown', function (e) {
      if (e.key === 'Enter') {
        var q = input.value.trim();
        if (q) location.href = withQsFallback('loci.html?q=' + encodeURIComponent(q));
      }
      if (e.key === 'Escape') { results.hidden = true; input.blur(); }
    });
    document.addEventListener('click', function (e) {
      if (!input.closest('.search-wrap').contains(e.target)) {
        results.hidden = true;
      }
    });
  }

  // ── Signal profile (SVG bar chart, 16 MAGMA p-values) ────

  function memberProfileSVG(gene) {
    var mp = DATA.memberProfiles && DATA.memberProfiles.genes && DATA.memberProfiles.genes[gene];
    if (!mp) return '';
    var signals = mp.signals || {};
    var underflow = new Set(mp.underflow || []);
    var vals = TRAIT_KEYS.map(function (t) {
      var p  = signals[t] != null ? signals[t] : 1;
      var lp = (p > 0 && p <= 1) ? -Math.log10(p) : 0;
      return { trait: t, p: p, lp: Math.min(15, lp), u: underflow.has(t) };
    });
    var maxLp = Math.ceil(Math.max(4, Math.max.apply(null, vals.map(function (v) { return v.lp; })) + 0.3));
    var mL = 28, mR = 6, mT = 14, mB = 32, barW = 13, barGap = 4, n16 = TRAIT_KEYS.length;
    var chartW = n16 * barW + (n16 - 1) * barGap, chartH = 88;
    var svgW = mL + chartW + mR, svgH = mT + chartH + mB, cb = mT + chartH;
    var s = [];
    s.push('<svg class="signal-svg" viewBox="0 0 ' + svgW + ' ' + svgH + '" width="100%" style="max-width:' + (svgW + 40) + 'px"' +
           ' aria-label="' + esc(gene) + ' MAGMA signal profile across 16 cardiovascular traits" role="img">');
    s.push('<line x1="' + mL + '" y1="' + mT + '" x2="' + mL + '" y2="' + cb + '" class="signal-axis"/>');
    var y01 = cb - (2 / maxLp) * chartH, y001 = cb - (3 / maxLp) * chartH;
    s.push('<line x1="' + mL + '" y1="' + y01.toFixed(1) + '" x2="' + (mL + chartW) + '" y2="' + y01.toFixed(1) + '" class="signal-t01"/>');
    s.push('<line x1="' + mL + '" y1="' + y001.toFixed(1) + '" x2="' + (mL + chartW) + '" y2="' + y001.toFixed(1) + '" class="signal-t001"/>');
    s.push('<text x="' + (mL - 3) + '" y="' + (cb + 4) + '" class="signal-ylab" text-anchor="end">0</text>');
    s.push('<text x="' + (mL - 3) + '" y="' + (y01 + 4).toFixed(1) + '" class="signal-ylab" text-anchor="end">2</text>');
    s.push('<text x="' + (mL - 3) + '" y="' + (y001 + 4).toFixed(1) + '" class="signal-ylab" text-anchor="end">3</text>');
    var midY = (mT + cb) / 2;
    s.push('<text x="9" y="' + midY.toFixed(1) + '" class="signal-ytitle" transform="rotate(-90 9 ' + midY.toFixed(1) + ')">&#8722;log&#8321;&#8320;p</text>');
    vals.forEach(function (v, i) {
      var bx = mL + i * (barW + barGap), bh = Math.max(0, (v.lp / maxLp) * chartH), by = cb - bh, rx = bh >= 3 ? 2 : 0;
      var cls = v.p < 0.001 ? 'signal-hi' : v.p < 0.01 ? 'signal-md' : 'signal-ns';
      s.push('<rect x="' + bx + '" y="' + (bh > 0 ? by : cb - 1).toFixed(1) + '" width="' + barW + '" height="' + Math.max(1, bh).toFixed(1) +
             '" rx="' + rx + '" class="' + cls + '" style="cursor:default" data-sig-trait="' + esc(traitLabel(v.trait)) + '" data-sig-p="' + v.p +
             '" data-sig-underflow="' + (v.u ? '1' : '0') + '">');
      s.push('<title>' + esc(traitLabel(v.trait)) + ': ' + (v.u ? 'P ≤ 1×10⁻³⁰⁰' : 'p = ' + v.p.toExponential(2)) + '</title>');
      s.push('</rect>');
      var lx = bx + barW / 2, ly = cb + 5;
      s.push('<text x="' + lx.toFixed(1) + '" y="' + ly + '" class="signal-tlab" text-anchor="end"' +
             ' transform="rotate(-45 ' + lx.toFixed(1) + ' ' + ly + ')">' + esc(TRAIT_SHORT[v.trait] || v.trait) + '</text>');
    });
    s.push('</svg>');
    return s.join('');
  }

  function renderMemberCompare(locus) {
    var mp = DATA.memberProfiles;
    if (!mp || !mp.by_locus || !mp.genes) return '';
    var members = mp.by_locus[locus.locus_id];
    if (!members || members.length < 2) return '';
    function pidFor(g) { return 'mprof-' + locus.locus_id + '-' + g.replace(/[^A-Za-z0-9_]/g, ''); }
    var rows = members.map(function (g) {
      var d = mp.genes[g];
      if (!d) return '';
      var repBadge = d.is_representative ? ' <span class="member-rep-badge" title="Statistical locus representative">rep</span>' : '';
      var pid = pidFor(g);
      var qtxt = (d.q_bh != null) ? Number(d.q_bh).toExponential(1) : '—';
      return [
        '<tr class="member-row">',
        '  <td class="member-gene"><span class="member-gene-sym">' + esc(g) + '</span>' + repBadge + '</td>',
        '  <td class="member-num">' + d.ctas.toFixed(2) + '</td>',
        '  <td class="member-num">' + qtxt + '</td>',
        '  <td class="member-num">' + d.m.toFixed(2) + '</td>',
        '  <td class="member-num">' + d.pr.toFixed(2) + '</td>',
        '  <td class="member-act"><button type="button" class="member-view-btn" data-mprof="' + pid + '"' +
             ' aria-expanded="false" aria-controls="' + pid + '">View profile</button></td>',
        '</tr>',
      ].join('\n');
    }).join('');
    // Profiles live OUTSIDE the (nowrap) table as full-width blocks, so the wide
    // SVG scales to the drawer and never forces body-level horizontal overflow.
    var profiles = members.map(function (g) {
      if (!mp.genes[g]) return '';
      var pid = pidFor(g);
      return '<div class="member-profile" id="' + pid + '" hidden>' +
             '<div class="member-profile-cap">' + esc(g) + ' raw 16-trait MAGMA profile</div>' +
             memberProfileSVG(g) + '</div>';
    }).join('');
    // Display-only "multiple signals" note for within-locus profile hotspots
    // (see data/loci.json multi_signal, sourced from the pleiotropic-hotspots
    // provenance CSV). Honest framing: candidate, not resolved — formal
    // independence would need fine-mapping, which is out of scope here.
    var multiNote = locus.multi_signal
      ? '  <p class="member-multi-note">Multiple signals — the member genes here fall into ≥2 distinct cross-trait profiles (candidate independent signals; compare the profiles below).</p>'
      : '';
    return [
      '<details class="member-compare">',
      '  <summary>Compare member-gene evidence (' + members.length + ')</summary>',
      multiNote,
      '  <div class="member-compare-scroll">',
      '  <table class="member-table" aria-label="Member-gene evidence for this locus">',
      '    <thead><tr><th scope="col">Gene</th><th scope="col">CTAS</th><th scope="col">q</th>' +
      '<th scope="col">M</th><th scope="col">PR</th><th scope="col"><span class="sr-only">Profile</span></th></tr></thead>',
      '    <tbody>' + rows + '</tbody>',
      '  </table>',
      '  </div>',
      '  <div class="member-profiles">' + profiles + '</div>',
      '</details>',
    ].join('\n');
  }

  function renderSignalProfile(locus) {
    var signals = DATA.traitSignals && DATA.traitSignals[locus.lead_gene];
    if (!signals) return '';

    // Traits whose MAGMA P underflowed to exact zero (floored to 1e-300 for plotting).
    // These display as a BOUND (P ≤ 1×10⁻³⁰⁰), never as an exact measured value.
    var underflow = new Set(signals['__underflow'] || []);

    // Build array of {trait, p, lp, u} in fixed canonical Table 1 order (no
    // significance re-sort), so every locus profile is directly comparable.
    var vals = TRAIT_KEYS.map(function (t) {
      var p  = signals[t] != null ? signals[t] : 1;
      var lp = (p > 0 && p <= 1) ? -Math.log10(p) : 0;
      return { trait: t, p: p, lp: Math.min(15, lp), u: underflow.has(t) };
    });

    // Y-axis ceiling: at least 4 (covers both threshold lines with margin)
    var maxLp = Math.ceil(Math.max(4, Math.max.apply(null, vals.map(function (v) { return v.lp; })) + 0.3));

    // Collect coloc-supported trait keys (PP4 ≥ 0.5 only; matches panel threshold)
    var colocTraits = new Set();
    if (locus.has_coloc) {
      (DATA.coloc || []).filter(function (c) { return c.locus_id === locus.locus_id; })
        .forEach(function (c) {
          (c.tissues || []).forEach(function (t) {
            if (t.max_pp4 >= 0.5 && t.best_trait) colocTraits.add(t.best_trait);
          });
        });
    }

    // SVG layout constants
    var mL = 28, mR = 6, mT = 14, mB = 32;
    var barW = 13, barGap = 4;
    var n16  = TRAIT_KEYS.length;
    var chartW = n16 * barW + (n16 - 1) * barGap;  // 268 px
    var chartH = 88;
    var svgW = mL + chartW + mR;   // 302 px
    var svgH = mT + chartH + mB;   // 134 px
    var cb = mT + chartH;          // chart-bottom Y

    var s = [];  // SVG markup parts

    s.push('<svg class="signal-svg" viewBox="0 0 ' + svgW + ' ' + svgH + '" width="100%"' +
           ' style="max-width:' + (svgW + 40) + 'px"' +
           ' aria-label="MAGMA signal profile across 16 cardiovascular traits" role="img">');

    // Y-axis spine
    s.push('<line x1="' + mL + '" y1="' + mT + '" x2="' + mL + '" y2="' + cb + '" class="signal-axis"/>');

    // Threshold lines at -log10(p)=2 and -log10(p)=3
    var y01  = cb - (2 / maxLp) * chartH;
    var y001 = cb - (3 / maxLp) * chartH;
    s.push('<line x1="' + mL + '" y1="' + y01.toFixed(1)  + '" x2="' + (mL+chartW) + '" y2="' + y01.toFixed(1)  + '" class="signal-t01"/>');
    s.push('<line x1="' + mL + '" y1="' + y001.toFixed(1) + '" x2="' + (mL+chartW) + '" y2="' + y001.toFixed(1) + '" class="signal-t001"/>');

    // Y-axis tick labels
    s.push('<text x="' + (mL-3) + '" y="' + (cb+4)             + '" class="signal-ylab" text-anchor="end">0</text>');
    s.push('<text x="' + (mL-3) + '" y="' + (y01+4).toFixed(1) + '" class="signal-ylab" text-anchor="end">2</text>');
    s.push('<text x="' + (mL-3) + '" y="' + (y001+4).toFixed(1)+ '" class="signal-ylab" text-anchor="end">3</text>');

    // Y-axis title (-log₁₀p), rotated
    var midY = (mT + cb) / 2;
    s.push('<text x="9" y="' + midY.toFixed(1) + '" class="signal-ytitle"' +
           ' transform="rotate(-90 9 ' + midY.toFixed(1) + ')">&#8722;log&#8321;&#8320;p</text>');

    // Bars, coloc dots, and trait labels
    vals.forEach(function (v, i) {
      var bx = mL + i * (barW + barGap);
      var bh = Math.max(0, (v.lp / maxLp) * chartH);
      var by = cb - bh;
      var rx = bh >= 3 ? 2 : 0;
      var cls = v.p < 0.001 ? 'signal-hi' : v.p < 0.01 ? 'signal-md' : 'signal-ns';

      // Bar — data attributes drive the custom floating tooltip;
      // <title> kept as accessible fallback for screen readers / keyboards.
      s.push('<rect x="' + bx + '" y="' + (bh > 0 ? by : cb - 1).toFixed(1) +
             '" width="' + barW + '" height="' + Math.max(1, bh).toFixed(1) +
             '" rx="' + rx + '" class="' + cls + '" style="cursor:default"' +
             ' data-sig-trait="' + esc(traitLabel(v.trait)) + '" data-sig-p="' + v.p + '"' +
             ' data-sig-underflow="' + (v.u ? '1' : '0') + '">');
      s.push('<title>' + esc(traitLabel(v.trait)) + ': ' + (v.u ? 'P ≤ 1×10⁻³⁰⁰' : 'p = ' + v.p.toExponential(2)) + '</title>');
      s.push('</rect>');

      // Coloc dot above bar
      if (colocTraits.has(v.trait)) {
        var dy = Math.max(mT + 5, by - 5);
        s.push('<circle cx="' + (bx + barW/2).toFixed(1) + '" cy="' + dy.toFixed(1) + '" r="2.5" class="signal-dot"/>');
      }

      // Trait label rotated -45° anchored at bar base
      var lx = bx + barW / 2;
      var ly = cb + 5;
      s.push('<text x="' + lx.toFixed(1) + '" y="' + ly + '" class="signal-tlab" text-anchor="end"' +
             ' transform="rotate(-45 ' + lx.toFixed(1) + ' ' + ly + ')">' +
             esc(TRAIT_SHORT[v.trait] || v.trait) + '</text>');
    });

    s.push('</svg>');

    // Legend
    var leg = [
      '<span class="signal-leg"><span class="signal-swatch signal-hi"></span>p&thinsp;&lt;&thinsp;0.001</span>',
      '<span class="signal-leg"><span class="signal-swatch signal-md"></span>p&thinsp;&lt;&thinsp;0.01</span>',
      '<span class="signal-leg"><span class="signal-swatch signal-ns"></span>ns</span>',
    ];
    if (colocTraits.size > 0) {
      leg.push('<span class="signal-leg">' +
               '<svg width="9" height="9" style="vertical-align:middle;margin-right:3px"><circle cx="4.5" cy="4.5" r="2.5" class="signal-dot"/></svg>' +
               'Coloc</span>');
    }

    return [
      '<div class="detail-section">',
      '  <div class="detail-section-title">' + esc(locus.lead_gene) + ' raw 16-trait MAGMA profile' + helpBtn('signal_profile', 'help-inline') + '</div>',
      '  <div class="signal-wrap">' + s.join('') + '</div>',
      '  <div class="signal-legend">' + leg.join('') + '</div>',
      '</div>',
    ].join('\n');
  }

  // ── Tier-grouped table rows ────────────────────────────────
  // Used when FS.activeFilter is set (route/domain landing) and no explicit
  // sort has been applied — renders tier-group header rows between tiers.

  function renderGroupedRows(loci) {
    // Sort by tier rank then ctas_score descending within each tier
    var sorted = loci.slice().sort(function (a, b) {
      var ta = TIER_ORDER[a.tier] != null ? TIER_ORDER[a.tier] : 99;
      var tb = TIER_ORDER[b.tier] != null ? TIER_ORDER[b.tier] : 99;
      if (ta !== tb) return ta - tb;
      return b.ctas_score - a.ctas_score;
    });

    var rows    = [];
    var curTier = null;

    sorted.forEach(function (locus) {
      if (locus.tier !== curTier) {
        curTier    = locus.tier;
        var tc = tierClass(locus.tier);
        var n  = sorted.filter(function (l) { return l.tier === locus.tier; }).length;
        rows.push(
          '<tr class="tier-group-hdr">' +
          '<td colspan="7"><div class="tier-group-inner">' +
          '<span class="tier ' + tc + '">' + esc(locus.tier) + '</span>' +
          '<span class="tier-group-n">' + n + ' loc' + (n === 1 ? 'us' : 'i') + '</span>' +
          '</div></td></tr>'
        );
      }
      rows.push(renderRow(locus));
    });

    return rows.join('');
  }

  // ── Loci Table Row ────────────────────────────────────────

  function renderRow(locus) {
    var tc    = tierClass(locus.tier);
    var qn    = tc === 'overlap' ? 'qt1' : tc === 'breadth' ? 'qt2' : 'qt3';
    var geneCount = (locus.n_genes === 1) ? '1 gene' : (locus.n_genes + ' genes');
    var pct   = scorePct(locus.ctas_score);
    var prior = priorLabel(locus);
    var pp3Entry = DATA.colocPP3Map[locus.locus_id];
    var coloc;
    if (locus.has_coloc) {
      var pp4Val = '<span class="coloc-pp4-val">PP4: ' + pp4Display(locus.locus_max_pp4) + '</span>';
      if (pp3Entry && pp3Entry.best_pp3 != null && pp3Entry.best_pp3 >= 0.5) {
        coloc = '<span class="coloc-vals">' + pp4Val +
                '<span class="coloc-sep">&middot;</span>' +
                '<span class="coloc-pp3-val">PP3: ' + pp3Entry.best_pp3.toFixed(3) + '</span></span>';
      } else {
        coloc = pp4Val;
      }
    } else if (pp3Entry && pp3Entry.coloc_class === 'PP3-only') {
      coloc = '<span class="coloc-pp3-val">PP3: ' + (pp3Entry.best_pp3 != null ? pp3Entry.best_pp3.toFixed(3) : '—') + '</span>';
    } else if (locus.coloc_class === 'not_assessed') {
      coloc = '<span class="coloc-no" title="Colocalisation was not assessed for this representative">Not assessed</span>';
    } else {
      coloc = '<span class="coloc-no">&#8212;</span>';
    }

    return [
      '<tr class="loci-row ' + qn + '" data-locus-id="' + locus.locus_id + '" tabindex="0" aria-controls="detail-panel" aria-expanded="false">',
      '  <td><div class="cell-gene">' + esc(locus.lead_gene) + '</div>' +
             '<div class="gene-count-sub">' + geneCount + '</div></td>',
      '  <td><span class="tier ' + tc + '">' + esc(locus.tier) + '</span></td>',
      '  <td>',
      '    <div class="score-cell">',
      '      <div class="score-bar-bg"><div class="score-bar-fill" style="width:' + pct + '%"></div></div>',
      '      <span class="score-value">' + locus.ctas_score.toFixed(1) + '</span>',
      '    </div>',
      '  </td>',
      '  <td><span class="count-pill">' + (locus.pr_breadth != null ? locus.pr_breadth.toFixed(2) : '—') + '</span></td>',
      '  <td><span class="family-cell' + (locus.family_raw_broad ? ' family-cell-broad' : '') + '">' + esc(familyRawDisplay(locus)) + '</span></td>',
      '  <td><span class="support support-' + prior.cls + '">' + esc(prior.text) + '</span></td>',
      '  <td>' + coloc + '</td>',
      '</tr>',
    ].join('\n');
  }

  function renderLociTable(loci) {
    var tbody = document.getElementById('loci-tbody');
    var count = document.getElementById('loci-count');
    if (!tbody) return;

    if (!loci.length) {
      tbody.innerHTML = '<tr><td colspan="7" class="table-empty">No loci match the current filters.</td></tr>';
      if (count) count.textContent = '0 loci';
      return;
    }

    // Grouped view when a domain/route filter is active and no explicit user sort
    if (FS.activeFilter && !SORT.col) {
      tbody.innerHTML = renderGroupedRows(loci);
    } else {
      tbody.innerHTML = loci.map(renderRow).join('');
    }
    if (count) {
      var total = DATA.loci.length;
      count.textContent = loci.length === total
        ? total + ' loci'
        : loci.length + ' of ' + total + ' loci';
    }
  }

  // ── Detail Panel ──────────────────────────────────────────

  // opts.deepLink = true  → called on fresh page load via ?locus= param;
  //                          use an instant single-scroll directly to the panel,
  //                          offset by the sticky nav height so the panel is fully
  //                          visible and not hidden behind the nav bar.
  // opts.deepLink = false → called by a row click while already on the page;
  //                          smooth-scroll the selected row into view, then the panel.
  function openDetail(locusId, opts) {
    opts = opts || {};
    var locus = DATA.loci.find(function (l) { return l.locus_id === locusId; });
    if (!locus) return;
    selectedLocusId = locusId;

    var selectedRow = null;
    document.querySelectorAll('.loci-row').forEach(function (r) {
      var match = r.dataset.locusId === locusId;
      r.classList.toggle('selected', match);
      r.setAttribute('aria-expanded', String(match));
      if (match) selectedRow = r;
    });

    renderDetailPanel(locus);
    var panel = document.getElementById('detail-panel');
    if (!panel) return;
    panel.classList.add('open');

    if (opts.deepLink) {
      // panel.classList.add('open') set display:block synchronously above, so
      // scrollIntoView can compute the panel's position right now.
      // scroll-margin-top on .detail-panel (in styles.css) offsets the sticky
      // nav so the panel header is not hidden behind it — no scrollBy needed.
      panel.scrollIntoView({ behavior: 'auto', block: 'start' });
    } else {
      // User row-click: smooth scroll — first the selected row into view (so the
      // user sees which locus was clicked), then the panel below it.
      setTimeout(function () {
        if (selectedRow) selectedRow.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
        setTimeout(function () {
          panel.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
        }, 180);
      }, 50);
    }
  }

  function closeDetail() {
    selectedLocusId = null;
    document.querySelectorAll('.loci-row').forEach(function (r) {
      r.classList.remove('selected');
      r.setAttribute('aria-expanded', 'false');
    });
    var panel = document.getElementById('detail-panel');
    if (panel) panel.classList.remove('open');
  }

  var TIER_COLORS = {
    'q<0.001': '#4A7A58',
    'q<0.01':  '#4A7FA5',
    'q<0.05':  '#7E8A9A',
  };

  // Numeric rank for q-tier sorting: lower = more stringent (q<0.001 best)
  var TIER_ORDER = { 'q<0.001': 0, 'q<0.01': 1, 'q<0.05': 2 };

  // Loci whose annotation notes draw on biological/literature context for traits
  // that are not in the locus's own CTAS-significant trait set (trait_list_p01).
  // This is expected — CTAS is a gene-level recurrence framework and may not
  // capture every trait with prior locus-level or regulatory support.
  // Where flagged, a light indicator links to the Method page explanation.
  // Broader-context annotation flags are a curated-pass item; empty for now.
  var BROADER_CONTEXT_IDS = new Set();

  // Compact drawer section: FinnGen external support.
  // Uses existing coloc-tissue-row styling; adds no table column, route, filter or badge.
  function renderEvidenceSection(locus) {
    // ── FinnGen ──
    var fgVal;
    if (locus.finngen_available === false) {
      fgVal = 'Not evaluable';
    } else {
      var n = locus.finngen_n_traits_p01 || 0;
      var nf = locus.finngen_n_families_p01 || 0;
      fgVal = nf + ' of 6 clinical families';
      if (locus.finngen_families) {
        fgVal += ' · ' + esc(String(locus.finngen_families).split(';').join(', '));
      }
      fgVal += '<br><span class="detail-evi-sub">' + n + ' of 10 endpoints at P&lt;0.01';
      if (locus.finngen_sig_traits) {
        var codes  = String(locus.finngen_sig_traits).split(';');
        var shorts = codes.map(function (c) { return (FINNGEN_LABELS[c] || [c])[0]; });
        var fulls  = codes.map(function (c) { return (FINNGEN_LABELS[c] || [c, c])[1]; });
        fgVal += ' · <span title="' + esc(fulls.join('; ')) + '">' + esc(shorts.join(', ')) + '</span>';
      }
      fgVal += '</span>';
    }

    function row(label, value, helpKey) {
      return '<div class="detail-coloc-tissue-row detail-evi-row">' +
             '<span class="detail-coloc-tissue-name">' + label +
               (helpKey ? helpBtn(helpKey, 'help-inline') : '') + '</span>' +
             '<span class="detail-evi-value">' + value + '</span></div>';
    }

    // ── Comparator studies: four individual studies plus their union ──
    var studies = [
      locus.in_qiao_2025  ? 'Qiao 2025'  : '',
      locus.in_liu_2023   ? 'Liu 2023'   : '',
      locus.in_foley_2021 ? 'Foley 2021' : '',
      locus.in_pleio      ? 'PLEIO'      : '',
    ].filter(Boolean);
    var cmpVal = studies.length
      ? (locus.n_prior_studies || studies.length) + ' of 4 studies · ' + esc(studies.join(', '))
      : 'None of 4 studies';

    return [
      '<div class="detail-section">',
      '  <div class="detail-section-title">External support</div>',
      '  <div class="detail-coloc-tissues">',
      row('Representative-gene FinnGen support', fgVal, 'finngen'),
      row('Comparator studies', cmpVal, 'prior'),
      '  </div>',
      '  <p class="detail-coloc-note">FinnGen provides external support for the statistical locus ' +
      'representative only; CTAS was not rerun or recalculated in FinnGen.</p>',
      '</div>',
    ].join('\n');
  }

  function renderDetailPanel(locus) {
    var panel = document.getElementById('detail-panel');
    if (!panel) return;

    // Tint the panel top border to match the locus tier
    panel.style.borderTopColor = TIER_COLORS[locus.tier] || 'var(--accent-lt)';

    var tc    = tierClass(locus.tier);
    var kw    = kwFor(locus);

    // Tier badge — clickable
    var tierRoute = TIER_ROUTE[locus.tier];
    var tierBadge = tierRoute
      ? '<a class="tier ' + tc + ' tier-lg" href="' + esc(withQsFallback('loci.html?route=' + encodeURIComponent(tierRoute))) + '">' + esc(locus.tier) + '</a>'
      : '<span class="tier ' + tc + ' tier-lg">' + esc(locus.tier) + '</span>';

    // Property badges
    var pp3Entry = DATA.colocPP3Map[locus.locus_id];
    var badges = [
      locus.ultra_stable     ? badgeLink('Ultra-stable', 'badge-stable', 'loci.html?route=ultra_stable')   : '',
      locus.is_canonical_cvd ? badge('Canonical CVD', 'badge-canon')                                       : '',
      locus.has_coloc        ? badgeLink('Coloc',       'badge-coloc',  'loci.html?route=coloc_supported') : '',
      (!locus.has_coloc && pp3Entry && pp3Entry.coloc_class === 'PP3-only') ? badge('Regional', 'badge-regional') : '',
      locus.not_reported     ? badgeLink('Not in 4 comparators', 'badge-novel', 'loci.html?route=not_reported') : '',
      locus.multi_signal     ? '<span class="badge badge-multi" title="Member genes fall into ≥2 distinct cross-trait profiles (candidate independent signals — open Compare member-gene evidence)">Multiple signals</span>' : '',
    ].filter(Boolean).join(' ');

    // Gene chips
    var geneChips = locus.genes.map(function (g) {
      var isLead = g === locus.lead_gene;
      return chip(g, 'chip-gene' + (isLead ? ' chip-gene-lead' : ''));
    }).join('');

    // Always-visible locked CTAS membership block. The complete set of
    // CTAS-significant member genes is essential data, shown by default near the
    // top of the drawer (not hidden inside optional Locus details). Any relevant
    // nearby/adjacent context gene is named once in the biological annotation prose,
    // so no separate context-gene line is rendered here.
    var memberGenesHtml = [
      '  <div class="detail-section detail-members">',
      '    <div class="detail-section-title">CTAS-significant genes in this locus (' + locus.n_genes + ')' +
             helpBtn('representative', 'help-inline') + '</div>',
      '    <div class="detail-chips detail-member-chips">' + geneChips + '</div>',
      renderMemberCompare(locus),
      '  </div>',
    ].join('\n');

    // Family chips — navigate via ?family= so the filter chips area shows
    // "Family: Coronary" rather than the raw search query
    var famChips = (locus.family_list_p01 || []).map(function (f) {
      return chipLink(familyLabel(f), 'chip-family', 'loci.html?family=' + encodeURIComponent(f));
    }).join('');

    // Trait chips — navigate via ?trait=
    var traitChips = (locus.trait_list_p01 || []).map(function (t) {
      return chipLink(traitLabel(t), 'chip-trait', 'loci.html?trait=' + encodeURIComponent(t));
    }).join('');

    // Prior study badges
    var studyBadges = [
      locus.in_qiao_2025  ? badge('Qiao 2025',  'badge-prior') : '',
      locus.in_liu_2023   ? badge('Liu 2023',   'badge-prior') : '',
      locus.in_foley_2021 ? badge('Foley 2021', 'badge-prior') : '',
      locus.in_pleio      ? badge('PLEIO',      'badge-prior') : '',
    ].filter(Boolean).join(' ');

    // Coloc section
    var colocHtml = '';
    if (locus.has_coloc) {
      var colocGeneChips = (locus.coloc_genes || []).map(function (g) {
        return chip(g, 'chip-gene');
      }).join('');

      var colocEntries = (DATA.coloc || []).filter(function (c) {
        return c.locus_id === locus.locus_id;
      });
      var tissueMap = {};
      colocEntries.forEach(function (entry) {
        (entry.tissues || []).forEach(function (t) {
          if (!tissueMap[t.tissue] || tissueMap[t.tissue].max_pp4 < t.max_pp4) {
            tissueMap[t.tissue] = { tissue: t.tissue, max_pp4: t.max_pp4, best_trait: entry.best_trait, tier: t.tier };
          }
        });
      });
      var tissueList = Object.values(tissueMap).sort(function (a, b) { return b.max_pp4 - a.max_pp4; });

      var tissueRowsHtml = '';
      if (tissueList.length > 0) {
        var rows = tissueList.map(function (t) {
          var tierCls = t.tier === 'high' ? 'detail-coloc-tier-high' : 'detail-coloc-tier-mod';
          return [
            '<div class="detail-coloc-tissue-row">',
            '<span class="detail-coloc-tissue-name">' + esc(tissueLabel(t.tissue)) + '</span>',
            '<span class="detail-coloc-tissue-trait">' + esc(traitLabel(t.best_trait || '—')) + '</span>',
            '<span class="detail-coloc-tissue-pp4 ' + tierCls + '"><span class="pp4-label">PP4</span> ' + pp4Display(t.max_pp4) + '</span>',
            '</div>',
          ].join('');
        }).join('\n');
        tissueRowsHtml = '<div class="detail-coloc-tissues">' + rows + '</div>';
      }

      // PP4 block — 3-stat grid (Max PP4 / Best tissue / Best trait), symmetric with PP3 block
      var pp4BlockHtml = [
        '<div class="detail-coloc-block">',
        '  <div class="detail-coloc-pp-type detail-coloc-pp4-type">Shared signal &middot; PP4 &#8805; 0.5</div>',
        '  <div class="detail-coloc-grid">',
        '    <div class="detail-stat"><span class="detail-stat-label">Max PP4</span>',
        '      <span class="detail-stat-value">' + pp4Display(locus.locus_max_pp4) + '</span></div>',
        '    <div class="detail-stat"><span class="detail-stat-label">Best tissue</span>',
        '      <span class="detail-stat-value">' + esc(tissueLabel(locus.locus_best_tissue || '—')) + '</span></div>',
        '    <div class="detail-stat"><span class="detail-stat-label">Best trait</span>',
        '      <span class="detail-stat-value">' + esc(traitLabel(locus.locus_best_trait || '—')) + '</span></div>',
        '  </div>',
        tissueRowsHtml,
        '</div>',
      ].join('\n');

      // PP3 block — shown if locus also has a meaningful PP3 signal
      var pp3BlockHtml = '';
      if (pp3Entry && pp3Entry.best_pp3 != null && pp3Entry.best_pp3 >= 0.5) {
        pp3BlockHtml = [
          '<div class="detail-coloc-block">',
          '  <div class="detail-coloc-pp-type detail-coloc-pp3-type">Regional overlap &middot; PP3 &#8805; 0.5</div>',
          '  <div class="detail-coloc-grid">',
          '    <div class="detail-stat"><span class="detail-stat-label">Max PP3</span>',
          '      <span class="detail-stat-value">' + pp3Entry.best_pp3.toFixed(3) + '</span></div>',
          '    <div class="detail-stat"><span class="detail-stat-label">Best tissue</span>',
          '      <span class="detail-stat-value">' + esc(tissueLabel(pp3Entry.best_pp3_tissue || '—')) + '</span></div>',
          '    <div class="detail-stat"><span class="detail-stat-label">Best trait</span>',
          '      <span class="detail-stat-value">' + esc(traitLabel(pp3Entry.best_pp3_trait || '—')) + '</span></div>',
          '  </div>',
          '</div>',
        ].join('\n');
      }

      // Explanatory note — wording depends on whether PP3 is also shown
      var colocNote = pp3BlockHtml
        ? 'PP4 &#8805; 0.5 supports a shared causal-variant model for the GWAS and eQTL signals; ' +
          'PP3 indicates both are present in the region but are more consistent with distinct local variants. ' +
          'Colocalisation support is not proof of a causal variant, causal gene or mechanism. ' +
          'Tissue rows show the context reaching PP4 &#8805; 0.5, with the top-scoring GWAS trait.'
        : 'PP4 &#8805; 0.5 supports a shared causal-variant model for the GWAS and eQTL signals; ' +
          'it is not proof of a causal variant, causal gene or mechanism. ' +
          'The row shows the tissue context reaching this threshold with its top-scoring GWAS trait. ' +
          'The full colocalisation scan covered more tissue&#8211;trait combinations.';

      // Coloc gene(s) — rendered once at section level, shown whenever any coloc data exists
      var colocGeneRowHtml = colocGeneChips
        ? '<div class="detail-coloc-gene-summary">' +
          '<span class="detail-coloc-gene-label">' +
          ((locus.coloc_genes || []).length === 1 ? 'Coloc gene' : 'Coloc genes') +
          '</span>' +
          colocGeneChips +
          '</div>'
        : '';

      colocHtml = [
        '<div class="detail-section">',
        '  <div class="detail-section-title">Colocalisation' + helpBtn('coloc_section', 'help-inline') + '</div>',
        pp4BlockHtml,
        pp3BlockHtml,
        colocGeneRowHtml,
        '  <p class="detail-coloc-note">' + colocNote + '</p>',
        '</div>',
      ].join('\n');
    } else if (pp3Entry && pp3Entry.coloc_class === 'PP3-only') {
      var pp3GeneChips = pp3Entry.best_pp3_gene
        ? chip(pp3Entry.best_pp3_gene, 'chip-gene')
        : '';
      var pp3GeneRowHtml = pp3GeneChips
        ? '<div class="detail-coloc-gene-summary">' +
          '<span class="detail-coloc-gene-label">Coloc gene</span>' +
          pp3GeneChips +
          '</div>'
        : '';
      colocHtml = [
        '<div class="detail-section">',
        '  <div class="detail-section-title">Colocalisation' + helpBtn('coloc_section', 'help-inline') + '</div>',
        '  <div class="detail-coloc-block">',
        '    <div class="detail-coloc-pp-type detail-coloc-pp3-type">Regional overlap &middot; PP3 &#8805; 0.5</div>',
        '    <div class="detail-coloc-grid">',
        '      <div class="detail-stat"><span class="detail-stat-label">Max PP3</span>',
        '        <span class="detail-stat-value">' + (pp3Entry.best_pp3 != null ? pp3Entry.best_pp3.toFixed(3) : '—') + '</span></div>',
        '      <div class="detail-stat"><span class="detail-stat-label">Best tissue</span>',
        '        <span class="detail-stat-value">' + esc(tissueLabel(pp3Entry.best_pp3_tissue || '—')) + '</span></div>',
        '      <div class="detail-stat"><span class="detail-stat-label">Best trait</span>',
        '        <span class="detail-stat-value">' + esc(traitLabel(pp3Entry.best_pp3_trait || '—')) + '</span></div>',
        '    </div>',
        '  </div>',
        pp3GeneRowHtml,
        '  <p class="detail-coloc-note">PP4 &#8805; 0.5 would support a shared causal-variant model for the GWAS and eQTL signals; ' +
           'PP3 indicates both signals are present in the region but are more consistent with distinct local variants. ' +
           'Colocalisation support is not proof of a causal variant, causal gene or mechanism. ' +
           'The values shown are the single strongest PP3 result from the full colocalisation scan.</p>',
        '</div>',
      ].join('\n');
    } else if (locus.coloc_class === 'not_assessed') {
      // Not assessed is missing evidence, never a negative result or PP4 = 0.
      var naReason = locus.coloc_disposition === 'not_assessed_MHC'
        ? 'the extended MHC region was excluded from colocalisation because its long-range linkage disequilibrium prevents interpretable single-variant tests'
        : 'no usable eQTL data were available for this representative gene in the tissues tested';
      colocHtml = [
        '<div class="detail-section">',
        '  <div class="detail-section-title">Colocalisation' + helpBtn('coloc_section', 'help-inline') + '</div>',
        '  <p class="detail-coloc-note">Not assessed: ' + naReason + '. This is missing evidence, not a negative colocalisation result.</p>',
        '</div>',
      ].join('\n');
    }

    // Annotation
    var annotHtml = '';
    if (locus.annotation_note) {
      // Evidence-context taxonomy (labels the locus-level evidence, not a CTAS tier
      // and not a causal-gene class). Legacy keys retained as a fallback.
      var EVIDENCE_LABELS = {
        'established_cross_cvd': 'Established cross-CVD',
        'established_cvd':       'Established CVD',
        'cvd_relevant_biology':  'CVD-relevant biology',
        'emerging_cross_trait':  'Emerging cross-trait',
        'known_cvd':              'Established CVD',
        'known_pleiotropic_cvd':  'Established cross-CVD',
        'cvd_adjacent_plausible': 'CVD-relevant biology',
        'emerging_ctas_specific': 'Emerging cross-trait',
      };
      var catKey = locus.evidence_category || locus.bio_category;
      var bioCat = catKey
        ? (EVIDENCE_LABELS[catKey] || catKey.replace(/_/g, ' '))
        : '';
      // Clickable, labelled sources kept behind an initially-collapsed "Sources (N)"
      // disclosure (PubMed/DOI for genetics/function; NCBI Gene/UniProt for basic
      // molecular function). Never raw unexplained URLs.
      var srcHtml = (locus.annotation_sources && locus.annotation_sources.length)
        ? '<details class="detail-disclosure detail-sources-disclosure">' +
          '<summary class="detail-disclosure-summary">Sources (' + locus.annotation_sources.length + ')</summary>' +
          '<div class="detail-disclosure-body detail-sources">' +
          locus.annotation_sources.map(function (s) {
            return '<a class="detail-source-link" href="' + esc(s.url) +
                   '" target="_blank" rel="noopener noreferrer">' + esc(s.label) + '</a>';
          }).join('<span class="detail-source-sep">·</span>') +
          '</div></details>'
        : '';
      // Light indicator for loci whose annotation note draws on broader
      // literature/biological context beyond the CTAS-significant trait set.
      var contextNote = BROADER_CONTEXT_IDS.has(locus.locus_id)
        ? '<div class="annotation-context-row">' +
            '<a href="method.html#annotation-context" class="annotation-context-link">' +
            'Annotation extends beyond the CTAS trait list' +
            '<svg width="10" height="10" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round" style="display:inline-block;vertical-align:middle;margin-left:3px;opacity:0.7"><path d="M18 13v6a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V8a2 2 0 0 1 2-2h6"/><polyline points="15 3 21 3 21 9"/><line x1="10" y1="14" x2="21" y2="3"/></svg>' +
            '</a></div>'
        : '';

      // The relevant biological-context gene is now named naturally, once, inside
      // annotation_note itself; the separate "Biological context gene:" line is no
      // longer rendered. The biological_context_gene field is retained in the data
      // (search keyword + build provenance) but not displayed as a distinct line.

      // Literature-derived prior context (curated known_family), shown with the biology,
      // omitted entirely where no defensible prior context exists (never "None").
      var priorCtxHtml = locus.known_family
        ? '  <p class="detail-priorctx"><span class="detail-priorctx-label">Prior cardiovascular context:</span> ' +
          '<span class="detail-priorctx-val">' + esc(familyLabel(locus.known_family)) + '</span>' +
          helpBtn('known_family', 'help-inline') + '</p>'
        : '';

      annotHtml = [
        '<div class="detail-section">',
        '  <div class="detail-section-title">Biological annotation</div>',
        '  <p class="detail-annot">' + esc(locus.annotation_note) + '</p>',
        priorCtxHtml,
        srcHtml,
        '  <div class="detail-chips" style="margin-top:8px;">',
        bioCat ? chip(bioCat, 'chip-bio') : '',
        studyBadges,
        '  </div>',
        contextNote,
        '</div>',
      ].join('\n');
    }

    // Position string — for single-gene point loci start === end, so
    // show just the coordinate rather than the degenerate "X–X" range.
    var posStr = (locus.start && locus.end && locus.start !== locus.end)
      ? (locus.start).toLocaleString() + '–' + (locus.end).toLocaleString()
      : (locus.start || 0).toLocaleString();

    panel.innerHTML = [
      '<div class="detail-panel-inner">',
      '  <button class="detail-close" aria-label="Close detail">&#215;</button>',

      '  <div class="detail-header">',
      '    <div>',
      '      <div class="detail-gene">' + esc(locus.lead_gene) + '</div>',
      // keyword chips removed from the header (keyword data stays searchable/filterable)
      '      <div class="detail-sub">' + esc(publicLocusLabel(locus.locus_id)) + ' &middot; Chr&nbsp;' + esc(locus.chr) +
             ':&nbsp;' + posStr + '</div>',
      '      <div class="detail-sub detail-rep-label">Statistical representative' + helpBtn('representative', 'help-inline') + '</div>',
      '    </div>',
      '    <div class="detail-badges">',
      '      ' + tierBadge,
      '      ' + badges,
      '    </div>',
      '  </div>',

      // Primary statistical summary: the CTAS decomposition (CTAS = M × log2(1+PR)),
      // the separate cross-family breadth (Family Neff / leading family), and stability.
      '  <div class="detail-section-title detail-rep-evidence-title">Representative-gene CTAS evidence' + helpBtn('rep_evidence', 'help-inline') + '</div>',
      '  <div class="detail-meta-grid">',
      '    <div class="detail-stat"><div class="detail-stat-head"><span class="detail-stat-label">CTAS</span>' + helpBtn('ctas') + '</div>',
      '      <span class="detail-stat-value">' + locus.ctas_score.toFixed(3) + '</span></div>',
      '    <div class="detail-stat"><div class="detail-stat-head"><span class="detail-stat-label">Magnitude, M</span>' + helpBtn('magnitude') + '</div>',
      '      <span class="detail-stat-value">' + (locus.magnitude_m != null ? locus.magnitude_m.toFixed(2) : '—') + '</span></div>',
      '    <div class="detail-stat"><div class="detail-stat-head"><span class="detail-stat-label">Effective breadth, PR</span>' + helpBtn('pr') + '</div>',
      '      <span class="detail-stat-value">' + (locus.pr_breadth != null ? locus.pr_breadth.toFixed(2) : '—') + '</span></div>',
      '    <div class="detail-stat"><div class="detail-stat-head"><span class="detail-stat-label">Family N&#8203;eff</span>' + helpBtn('family_neff') + '</div>',
      '      <span class="detail-stat-value">' + (locus.neff_family_breadth != null ? Number(locus.neff_family_breadth).toFixed(2) : '—') + '</span></div>',
      // Two calculated family summaries: raw-profile (primary) and dependence-corrected (secondary).
      // The literature-derived "prior cardiovascular context" (known_family) is shown separately
      // inside the Biological annotation section, not alongside these calculated summaries.
      '    <div class="detail-stat"><div class="detail-stat-head"><span class="detail-stat-label">Raw-profile leading family</span>' + helpBtn('family_raw') + '</div>',
      '      <span class="detail-stat-value' + (locus.family_raw_broad ? ' detail-stat-value-muted' : '') + '">' + esc(familyRawDisplay(locus)) + '</span></div>',
      '    <div class="detail-stat detail-stat-secondary"><div class="detail-stat-head"><span class="detail-stat-label">Leading family after dependence correction</span>' + helpBtn('family_ctas_weighted') + '</div>',
      '      <span class="detail-stat-value">' + esc(locus.leading_family_display || familyLabel(locus.largest_share_family)) + '</span>' +
             (locus.leading_family_pct_display ? '<span class="detail-stat-sub">' + esc(locus.leading_family_pct_display) + '</span>' : '') + '</div>',
      '    <div class="detail-stat"><div class="detail-stat-head"><span class="detail-stat-label">Stability</span>' + helpBtn('stability') + '</div>',
      '      <span class="detail-stat-value">' + esc(locus.stability_class || '—') +
             ' (' + (locus.loo_runs_retained || 0) + '/16)</span></div>',
      '  </div>',

      memberGenesHtml,

      renderSignalProfile(locus),

      // Descriptive context only: raw named-trait / family counts at nominal
      // thresholds. These are a compact summary of the profile above, NOT
      // components of CTAS and not a substitute for the effective breadth PR.
      '  <div class="detail-section">',
      '    <div class="detail-section-title">Raw named-trait summary <span class="detail-context-tag">descriptive context</span></div>',
      '    <div class="detail-meta-grid detail-meta-grid-secondary">',
      '      <div class="detail-stat"><div class="detail-stat-head"><span class="detail-stat-label">Traits P&lt;0.01</span>' + helpBtn('traits_p01') + '</div>',
      '        <span class="detail-stat-value">' + locus.n_traits_p01 + '</span></div>',
      '      <div class="detail-stat"><div class="detail-stat-head"><span class="detail-stat-label">Traits P&lt;0.001</span>' + helpBtn('traits_p001') + '</div>',
      '        <span class="detail-stat-value">' + locus.n_traits_p001 + '</span></div>',
      '      <div class="detail-stat"><div class="detail-stat-head"><span class="detail-stat-label">Families P&lt;0.01</span>' + helpBtn('families_p01') + '</div>',
      '        <span class="detail-stat-value">' + locus.n_families_p01 + '</span></div>',
      '      <div class="detail-stat"><div class="detail-stat-head"><span class="detail-stat-label">Span</span>' + helpBtn('span') + '</div>',
      '        <span class="detail-stat-value">' + formatSpan(locus.span_kb) + '</span></div>',
      '    </div>',
      // raw-trait prose removed: each box carries its own `?` definition (no duplication)
      '  </div>',

      annotHtml,

      // Progressive disclosure: detailed gene/family/named-trait lists (initially collapsed).
      // Keyword chips are no longer rendered here; the data remains searchable/filterable.
      '  <details class="detail-section detail-disclosure">',
      '    <summary class="detail-disclosure-summary">Locus details' +
             ' <span class="detail-disclosure-hint">families · traits</span></summary>',
      '    <div class="detail-disclosure-body">',
      '      <div class="detail-subsection-title">Families at P&lt;0.01</div>',
      '      <div class="detail-chips">' + famChips + '</div>',
      '      <div class="detail-subsection-title">Traits at P&lt;0.01</div>',
      '      <div class="detail-chips">' + traitChips + '</div>',
      '    </div>',
      '  </details>',

      colocHtml,

      renderEvidenceSection(locus),

      '  <div class="detail-method-link">',
      '    <a href="method.html">Method &amp; field definitions &#8594;</a>',
      '  </div>',

      '</div>',
    ].join('\n');

    panel.querySelector('.detail-close').addEventListener('click', closeDetail);

    // Member-gene "View profile" disclosure: toggle the target profile; keep at
    // most one open at a time within the drawer. Keyboard-operable (real button).
    // The panel element persists across renders (innerHTML is replaced), so this
    // delegated listener is attached ONCE — attaching it per render would stack
    // handlers and make each click toggle N times (net no-op on even counts,
    // which is why some drawers' profiles appeared to never open).
    if (!panel.dataset.memberHandlerAttached) {
      panel.dataset.memberHandlerAttached = '1';
      panel.addEventListener('click', function (e) {
        var btn = e.target.closest ? e.target.closest('.member-view-btn') : null;
        if (!btn) return;
        var id = btn.getAttribute('data-mprof');
        var target = document.getElementById(id);
        if (!target) return;
        var willOpen = target.hidden;
        panel.querySelectorAll('.member-profile').forEach(function (p) { p.hidden = true; });
        panel.querySelectorAll('.member-view-btn').forEach(function (b) {
          b.setAttribute('aria-expanded', 'false'); b.textContent = 'View profile';
        });
        if (willOpen) {
          target.hidden = false;
          btn.setAttribute('aria-expanded', 'true'); btn.textContent = 'Hide profile';
        }
      });
    }
  }

  // ── Active filter chips area ──────────────────────────────
  // Renders named, individually-removable chips for every active filter.
  // Called after every filter state change.
  function renderFilterChips() {
    var area = document.getElementById('filter-chips-area');
    var list = document.getElementById('filter-chips-list');
    if (!area || !list) return;

    var chips = [];

    // 1. Route / domain (primary context filter)
    if (FS.activeFilter) {
      var ctxType = FS.activeFilter.type === 'domain' ? 'Domain' : 'Route';
      chips.push({ type: ctxType, value: FS.activeFilter.label, key: 'active' });
    }

    // 2. Family / trait / text search
    if (FS.query) {
      if (FS.filterLabel) {
        chips.push({
          type:  FS.filterLabel.type === 'family' ? 'Family' : 'Trait',
          value: FS.filterLabel.label,
          key:   'query'
        });
      } else {
        chips.push({ type: 'Search', value: '“' + FS.query + '”', key: 'query' });
      }
    }

    // 3. q-threshold chip (nested; one at a time)
    if (FS.qThreshold) {
      var qLbl = { q05: 'q<0.05', q01: 'q<0.01', q001: 'q<0.001' }[FS.qThreshold];
      chips.push({ type: 'q threshold', value: qLbl, key: 'qthr' });
    }

    // 4. Coloc / comparator toggles
    if (FS.colocOnly) chips.push({ type: 'Coloc', value: 'Supported', key: 'coloc' });
    if (FS.novelOnly) chips.push({ type: 'Not in 4 comparators', value: 'only', key: 'novel' });

    if (!chips.length) { area.hidden = true; return; }

    area.hidden = false;
    list.innerHTML = chips.map(function (c) {
      return [
        '<span class="fac">',
        '<span class="fac-type">' + esc(c.type) + ':</span>',
        '<span class="fac-value">' + esc(c.value) + '</span>',
        '<button class="fac-remove" data-remove="' + esc(c.key) + '"',
        '        aria-label="Remove ' + esc(c.type) + ' filter">&#215;</button>',
        '</span>',
      ].join('');
    }).join('');
  }

  // ── Sort headers (loci page) ──────────────────────────────

  function initSort() {
    // Columns that default to descending on first click (numeric — "best first")
    var descFirst = new Set(['ctas_score', 'n_traits_p01', 'n_families_p01', 'prior', 'coloc']);
    // 'tier' sorts asc first (q<0.001 = rank 0, so low->high = most stringent first)
    // 'lead_gene' sorts asc first (A→Z)

    document.querySelectorAll('th[data-sort-col]').forEach(function (th) {
      var col = th.dataset.sortCol;
      th.setAttribute('tabindex', '0');
      th.setAttribute('role', 'columnheader');
      th.setAttribute('aria-sort', 'none');

      function doSort() {
        var newDir;
        if (SORT.col === col) {
          newDir = -SORT.dir;           // same column: flip direction
        } else {
          newDir = descFirst.has(col) ? -1 : 1;  // new column: use default direction
        }
        SORT.col = col;
        SORT.dir = newDir;

        // Update visual state on all headers
        document.querySelectorAll('th[data-sort-col]').forEach(function (h) {
          h.classList.remove('sort-asc', 'sort-desc');
          h.setAttribute('aria-sort', 'none');
        });
        th.classList.add(SORT.dir === 1 ? 'sort-asc' : 'sort-desc');
        th.setAttribute('aria-sort', SORT.dir === 1 ? 'ascending' : 'descending');

        applyFilters();
      }

      th.addEventListener('click', doSort);
      th.addEventListener('keydown', function (e) {
        if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); doSort(); }
      });
    });
  }

  // ── Filters (loci page) ───────────────────────────────────

  function filteredLoci() {
    return DATA.loci.filter(function (l) {
      if (FS.locusIds && !FS.locusIds.has(l.locus_id))  return false;
      if (FS.qThreshold) {                                // nested q-threshold (survival flags)
        var flag = FS.qThreshold === 'q001' ? 'survives_q001'
                 : FS.qThreshold === 'q01'  ? 'survives_q01'
                 :                            'survives_q05';
        if (!l[flag]) return false;
      }
      if (FS.colocOnly && !l.has_coloc)                   return false;
      if (FS.novelOnly && !l.not_reported)                return false;
      if (FS.query) {
        var q  = FS.query.toLowerCase();
        var kw = kwFor(l);
        var colocTerms = colocSearchTerms(l);
        return l.lead_gene.toLowerCase().includes(q)
          || (l.annotation_note && l.annotation_note.toLowerCase().includes(q))
          || l.genes.some(function (g) { return g.toLowerCase().includes(q); })
          || (l.family_list_p01 || []).some(function (f) {
               return f.toLowerCase().includes(q) || familyLabel(f).toLowerCase().includes(q);
             })
          || (l.leading_family_keys || []).some(function (k) {
               return k.toLowerCase().includes(q) || familyLabel(k).toLowerCase().includes(q);
             })
          || (l.trait_list_p01 || []).some(function (t) {
               return t.toLowerCase().includes(q) || traitLabel(t).toLowerCase().includes(q);
             })
          || (kw && (kw.search_keywords || []).some(function (k) { return k.toLowerCase().includes(q); }))
          || (kw && (kw.biological_context || []).some(function (k) { return k.toLowerCase().includes(q); }))
          || (kw && kw.locus_candidate && kw.locus_candidate.toLowerCase().includes(q))
          || colocTerms.some(function (t) { return String(t).toLowerCase().includes(q); })
          || evidenceSearchTerms(l).some(function (t) { return t.toLowerCase().includes(q); });
      }
      return true;
    });
  }

  // ── Sort ──────────────────────────────────────────────────
  // Returns a sorted copy of the loci array according to SORT state.
  // Called after filteredLoci() so sort is always applied on top of filters.
  function sortedLoci(loci) {
    if (!SORT.col || !SORT.dir) return loci;
    var col = SORT.col;
    var dir = SORT.dir;

    return loci.slice().sort(function (a, b) {
      var va, vb;
      switch (col) {
        case 'lead_gene':
          va = a.lead_gene; vb = b.lead_gene; break;
        case 'tier':
          // Sort by q-tier rank (q<0.001 = 0 < q<0.01 = 1 < q<0.05 = 2)
          va = TIER_ORDER[a.tier] != null ? TIER_ORDER[a.tier] : 99;
          vb = TIER_ORDER[b.tier] != null ? TIER_ORDER[b.tier] : 99;
          break;
        case 'ctas_score':
          va = a.ctas_score; vb = b.ctas_score; break;
        case 'pr_breadth':
          va = a.pr_breadth != null ? a.pr_breadth : -1;
          vb = b.pr_breadth != null ? b.pr_breadth : -1; break;
        case 'top_family':
          // sort by the raw primary family; "Broad" loci sort to the end
          va = a.family_raw_broad ? '~' : familyRawDisplay(a);
          vb = b.family_raw_broad ? '~' : familyRawDisplay(b); break;
        case 'prior':
          // n_prior_studies: 0 = not reported, 1–4 = comparator studies matched
          va = a.n_prior_studies || 0; vb = b.n_prior_studies || 0; break;
        case 'coloc':
          // Ordered categorical: PP4-supported (2+pp4) > PP3-only (1+pp3) > none (0) > untested (-1)
          va = a.has_coloc ? 2 + (a.locus_max_pp4 || 0)
             : (DATA.colocPP3Map[a.locus_id] && DATA.colocPP3Map[a.locus_id].coloc_class === 'PP3-only')
               ? 1 + (DATA.colocPP3Map[a.locus_id].best_pp3 || 0)
             : DATA.colocPP3Map[a.locus_id] ? 0 : -1;
          vb = b.has_coloc ? 2 + (b.locus_max_pp4 || 0)
             : (DATA.colocPP3Map[b.locus_id] && DATA.colocPP3Map[b.locus_id].coloc_class === 'PP3-only')
               ? 1 + (DATA.colocPP3Map[b.locus_id].best_pp3 || 0)
             : DATA.colocPP3Map[b.locus_id] ? 0 : -1;
          break;
        default: return 0;
      }
      if (typeof va === 'string') return dir * va.localeCompare(vb);
      return dir * (va - vb);
    });
  }

  function applyFilters() {
    renderLociTable(sortedLoci(filteredLoci()));
    closeDetail();
    renderFilterChips();
  }

  function initLociFilters() {
    var input    = document.getElementById('loci-filter');
    var colocBtn = document.getElementById('filter-coloc');
    var novelBtn = document.getElementById('filter-novel');

    // ── Text input ──
    if (input) {
      input.addEventListener('input', function () {
        FS.query       = input.value.trim();
        FS.filterLabel = null;   // manual typing clears any semantic label
        applyFilters();
      });
    }

    // ── q-threshold buttons (nested; one active at a time) ──
    document.querySelectorAll('[data-tier-filter]').forEach(function (btn) {
      btn.addEventListener('click', function () {
        var v = btn.dataset.tierFilter;
        if (FS.qThreshold === v) {
          FS.qThreshold = null;
          btn.classList.remove('active');
        } else {
          FS.qThreshold = v;
          document.querySelectorAll('[data-tier-filter]').forEach(function (b) {
            b.classList.toggle('active', b.dataset.tierFilter === v);
          });
        }
        applyFilters();
      });
    });

    // ── Coloc / Novel toggles ──
    if (colocBtn) {
      colocBtn.addEventListener('click', function () {
        FS.colocOnly = !FS.colocOnly;
        colocBtn.classList.toggle('active', FS.colocOnly);
        applyFilters();
      });
    }
    if (novelBtn) {
      novelBtn.addEventListener('click', function () {
        FS.novelOnly = !FS.novelOnly;
        novelBtn.classList.toggle('active', FS.novelOnly);
        applyFilters();
      });
    }

    // ── Row click ──
    var tbody = document.getElementById('loci-tbody');
    if (tbody) {
      tbody.addEventListener('click', function (e) {
        var row = e.target.closest('.loci-row');
        if (!row) return;
        var id = row.dataset.locusId;
        if (id === selectedLocusId) closeDetail();
        else openDetail(id);
      });
      // Keyboard: Enter / Space on a focused row toggles its detail drawer.
      tbody.addEventListener('keydown', function (e) {
        if (e.key !== 'Enter' && e.key !== ' ') return;
        var row = e.target.closest('.loci-row');
        if (!row) return;
        e.preventDefault();
        var id = row.dataset.locusId;
        if (id === selectedLocusId) closeDetail();
        else openDetail(id);
      });
    }

    // ── Filter chips — chip removal event delegation ──
    var chipsList = document.getElementById('filter-chips-list');
    if (chipsList) {
      chipsList.addEventListener('click', function (e) {
        var btn = e.target.closest('[data-remove]');
        if (!btn) return;
        var key = btn.dataset.remove;

        if (key === 'active') {
          FS.locusIds     = null;
          FS.activeFilter = null;
          if (!FS.query && !FS.qThreshold && !FS.colocOnly && !FS.novelOnly) {
            resetPageHeader();
          }
          var u = new URL(location.href);
          u.searchParams.delete('domain');
          u.searchParams.delete('route');
          clearQsHash(u); history.replaceState(null, '', u.toString());

        } else if (key === 'query') {
          FS.query       = '';
          FS.filterLabel = null;
          if (input) input.value = '';
          if (!FS.activeFilter) resetPageHeader();
          var u = new URL(location.href);
          u.searchParams.delete('q');
          u.searchParams.delete('family');
          u.searchParams.delete('trait');
          clearQsHash(u); history.replaceState(null, '', u.toString());

        } else if (key === 'qthr') {
          FS.qThreshold = null;
          document.querySelectorAll('[data-tier-filter]').forEach(function (b) {
            b.classList.remove('active');
          });

        } else if (key === 'coloc') {
          FS.colocOnly = false;
          if (colocBtn) colocBtn.classList.remove('active');

        } else if (key === 'novel') {
          FS.novelOnly = false;
          if (novelBtn) novelBtn.classList.remove('active');
        }

        applyFilters();
      });
    }

    // ── Clear-all button ──
    var clearAllBtn = document.getElementById('filter-clear-all');
    if (clearAllBtn) {
      clearAllBtn.addEventListener('click', function () {
        FS.query        = '';
        FS.filterLabel  = null;
        FS.qThreshold   = null;
        FS.colocOnly    = false;
        FS.novelOnly    = false;
        FS.locusIds     = null;
        FS.activeFilter = null;

        if (input) input.value = '';
        document.querySelectorAll('[data-tier-filter]').forEach(function (b) { b.classList.remove('active'); });
        if (colocBtn) colocBtn.classList.remove('active');
        if (novelBtn) novelBtn.classList.remove('active');

        resetPageHeader();

        var u = new URL(location.href);
        ['domain','route','q','family','trait','tier','coloc'].forEach(function(k){ u.searchParams.delete(k); });
        clearQsHash(u); history.replaceState(null, '', u.toString());

        applyFilters();
      });
    }

    // ── Process URL params ──
    var params = currentParams();

    // Generic text search (?q=)
    if (params.has('q') && input) {
      FS.query    = params.get('q');
      input.value = FS.query;
    }

    // q-threshold pre-filter (?tier=q05|q01|q001 — nested, one at a time)
    if (params.has('tier')) {
      var tv = params.get('tier').trim();
      if (['q05', 'q01', 'q001'].indexOf(tv) !== -1) {
        FS.qThreshold = tv;
        document.querySelectorAll('[data-tier-filter]').forEach(function (b) {
          b.classList.toggle('active', b.dataset.tierFilter === tv);
        });
      }
    }

    // Family chip nav (?family=Coronary) — semantic filter with display label
    if (params.has('family')) {
      var fk = params.get('family');
      FS.query       = fk;
      FS.filterLabel = { type: 'family', label: familyLabel(fk) };
      var fl         = familyLabel(fk);
      setNavTitle(fl);
      setPageHeader('Disease family', fl + ' loci', fl + ' loci among the 91 CTAS loci — click any row for details.');
    }

    // Trait chip nav (?trait=CAD_I25) — semantic filter with display label
    if (params.has('trait')) {
      var tk = params.get('trait');
      FS.query       = tk;
      FS.filterLabel = { type: 'trait', label: traitLabel(tk) };
      var tl         = traitLabel(tk);
      setNavTitle(tl);
      setPageHeader('Trait', tl + ' loci', 'Loci with signal at p < 0.01 for ' + tl + ' — click any row for details.');
    }

    // Domain pre-filter (?domain=)
    if (params.has('domain')) {
      var domainId = params.get('domain');
      var domain   = DATA.domains.find(function (d) { return d.id === domainId; });
      if (domain) {
        FS.locusIds     = new Set([].concat(domain.locus_ids));
        FS.activeFilter = { type: 'domain', label: domain.label };
        setNavTitle(domain.label);
        setPageHeader('Disease domain', domain.label + ' loci', domain.label + ' loci among the 91 CTAS loci — click any row for details.');
      }
    }

    // Route pre-filter (?route=)
    if (params.has('route')) {
      var routeId = params.get('route');
      var route   = DATA.featured.find(function (r) { return r.id === routeId; });
      if (route) {
        FS.locusIds     = new Set([].concat(route.locus_ids));
        FS.activeFilter = { type: 'route', label: route.label };
        setNavTitle(route.label);
        setPageHeader('Featured route', route.label, route.description);
      }
    }

    // Coloc pre-filter (?coloc=)
    if (params.has('coloc') && colocBtn) {
      FS.colocOnly = true;
      colocBtn.classList.add('active');
    }

    // Novel pre-filter (?novel=)
    if (params.has('novel') && novelBtn) {
      FS.novelOnly = true;
      novelBtn.classList.add('active');
    }

    applyFilters();

    // Auto-open detail panel.
    // Priority: an explicit ?locus=V2Lxxx deep link takes precedence; otherwise
    // the homepage landscape click handler stored the target id in sessionStorage
    // before navigating here. Both are read once. Locus ids are strings.
    var pendingLocus = null;
    if (params.has('locus')) {
      pendingLocus = params.get('locus');
      sessionStorage.removeItem('ctas_open_locus');   // param wins over any stored id
    } else {
      pendingLocus = sessionStorage.getItem('ctas_open_locus');
    }
    if (pendingLocus) {
      sessionStorage.removeItem('ctas_open_locus');
      var pendingId = pendingLocus;
      if (pendingId && DATA.loci.some(function (l) { return l.locus_id === pendingId; })) {
        // Defer one tick so the browser finishes its own initial paint before
        // we scroll — keeps the transition visually clean.
        setTimeout(function () {
          openDetail(pendingId, { deepLink: true });
        }, 0);
      }
    }
  }

  // ── Page init ─────────────────────────────────────────────

  function initIndex() {
    renderLandscapePlot(DATA.loci);
    renderFeaturedCards(DATA.featured);
    renderDomainCards(DATA.domains);
    initHeroSearch(DATA.loci);
  }

  // ── Signal-profile bar tooltip ────────────────────────────
  // Delegated listener on #detail-panel; reads data-sig-* attributes placed
  // on each <rect> by renderSignalProfile.  One persistent tooltip div is
  // created at init time and reused for every hover.

  function initSignalTooltip() {
    var panel = document.getElementById('detail-panel');
    if (!panel) return;

    var tipEl = document.createElement('div');
    tipEl.className = 'sig-tooltip';
    tipEl.hidden = true;
    document.body.appendChild(tipEl);

    panel.addEventListener('mousemove', function (e) {
      var bar = e.target;
      if (!bar || !bar.getAttribute || !bar.getAttribute('data-sig-trait')) {
        tipEl.hidden = true;
        return;
      }
      var traitName = bar.getAttribute('data-sig-trait');
      var p         = parseFloat(bar.getAttribute('data-sig-p'));
      var isUnder   = bar.getAttribute('data-sig-underflow') === '1';
      var lp        = (p > 0 && p <= 1) ? (-Math.log10(p)).toFixed(2) : '0.00';
      // Underflowed MAGMA P (exact zero, floored to 1e-300 for plotting) is shown as a
      // bound, not an exact value; the \u2212log\u2081\u2080p line likewise reads as a lower bound.
      var pLine  = isUnder ? 'P\u2009\u2264\u20091\u00d710\u207b\u00b3\u2070\u2070'
                           : 'p\u2009=\u2009' + p.toExponential(2);
      var lpLine = isUnder ? '\u2212log\u2081\u2080p\u2009\u2265\u2009300'
                           : '\u2212log\u2081\u2080p\u2009=\u2009' + lp;
      tipEl.innerHTML =
        '<div class="sig-tip-trait">' + esc(traitName) + '</div>' +
        '<div class="sig-tip-vals">' +
          '<span>' + pLine + '</span>' +
          '<span class="sig-tip-lp">' + lpLine + '</span>' +
        '</div>';
      var tx = e.clientX + 14;
      var ty = e.clientY - 14;
      if (tx + 190 > window.innerWidth)  tx = e.clientX - 190 - 10;
      if (ty + 72  > window.innerHeight) ty = e.clientY - 72;
      tipEl.style.left = tx + 'px';
      tipEl.style.top  = ty + 'px';
      tipEl.hidden = false;
    });

    panel.addEventListener('mouseleave', function () { tipEl.hidden = true; });

    // Hide when user clicks outside the panel.
    document.addEventListener('click', function (e) {
      if (!panel.contains(e.target)) tipEl.hidden = true;
    });
  }

  function initLoci() {
    SCORE_MAX = DATA.loci.reduce(function (m, l) { return Math.max(m, l.ctas_score); }, 1);
    initSort();
    initLociFilters();
    initSignalTooltip();
  }

  async function init() {
    // Prevent the browser's automatic scroll restoration from overriding
    // programmatic scrolls on deep-link page loads (e.g. ?locus= arrivals).
    if ('scrollRestoration' in history) history.scrollRestoration = 'manual';

    setActiveNav();
    initMobileNav();
    initBackButton();

    var page = location.pathname.split('/').pop().replace(/\.html$/, '') || 'index';

    try {
      if (page === '' || page === 'index') {
        var results = await Promise.all([
          loadJSON('data/domains.json'),
          loadJSON('data/featured.json'),
          loadJSON('data/loci.json'),
          loadJSON('data/locus_keywords.json').catch(function () { return {}; }),
        ]);
        DATA.domains  = results[0];
        DATA.featured = results[1];
        DATA.loci     = normalizeLoci(results[2]);
        DATA.keywords = results[3] || {};
        joinKeywords(DATA.loci, DATA.keywords);
        initIndex();

      } else if (page === 'loci') {
        var results = await Promise.all([
          loadJSON('data/loci.json'),
          loadJSON('data/domains.json'),
          loadJSON('data/featured.json'),
          loadJSON('data/coloc.json'),
          loadJSON('data/trait_signals.json').catch(function () { return {}; }),
          loadJSON('data/locus_keywords.json').catch(function () { return {}; }),
          loadJSON('data/coloc_pp3.json').catch(function () { return []; }),
          loadJSON('data/member_gene_profiles.json').catch(function () { return null; }),
        ]);
        DATA.loci         = normalizeLoci(results[0]);
        DATA.domains      = results[1];
        DATA.featured     = results[2];
        DATA.coloc        = results[3] || [];
        DATA.traitSignals = results[4] || {};
        DATA.keywords     = results[5] || {};
        DATA.colocPP3     = results[6] || [];
        DATA.memberProfiles = results[7] || null;
        DATA.colocPP3Map  = {};
        DATA.colocPP3.forEach(function (entry) {
          DATA.colocPP3Map[entry.locus_id] = entry;
        });
        joinKeywords(DATA.loci, DATA.keywords);
        initLoci();
      }

    } catch (err) {
      console.error('CTAS Explorer: data load failed —', err.message);
      var grid = document.getElementById('domain-grid') || document.getElementById('loci-tbody');
      if (grid) {
        grid.innerHTML = '<p style="padding:24px;color:var(--text-3)">Data could not be loaded. ' +
          'Serve the explorer via a local server (e.g. <code>npx serve .</code>).</p>';
      }
    }
  }

  document.addEventListener('DOMContentLoaded', init);

}());
