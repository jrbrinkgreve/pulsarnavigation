// Pulsar navigation docs: sidebar, prev/next links, theme toggle and equations.
// The page list lives here only, so adding a page = one line in PAGES.
// Works from file:// (no server, no external resources; KaTeX is in vendor/).

const PAGES = [
  { group: "Start" },
  { file: "index.html",       title: "Home" },
  { file: "physics.html",     title: "Physics background" },
  { file: "localisation.html", title: "Localisation: frames, time, relativity" },
  { file: "conventions.html", title: "Conventions and why" },
  { file: "pipeline.html",    title: "Pipeline overview" },
  { file: "parameters.html",  title: "Parameters" },
  { group: "Stages" },
  { file: "sky.html",         title: "1–2 Sky: pulsar and dispersion" },
  { file: "receiver.html",    title: "3–4 Receiver: noise, RFI, IQ, clock" },
  { file: "fullband.html",    title: "5–6 Full-band front end" },
  { file: "channels.html",    title: "5a–6b Channelized front end" },
  { file: "excision.html",    title: "RFI excision" },
  { file: "fold.html",        title: "7 Phase prediction and folding" },
  { file: "simulator.html",   title: "1–7 Fast simulator (long runs)" },
  { file: "toa.html",         title: "8 TOA estimation" },
  { file: "detection.html",   title: "Detection and good TOAs" },
  { file: "validation.html",  title: "9 Validation and checks" },
  { file: "navigation.html",  title: "10 Navigation" },
  { group: "Running" },
  { file: "scripts.html",     title: "Run scripts" },
  { file: "tests.html",       title: "Unit tests" },
  { group: "Results and plans" },
  { file: "history.html",     title: "Development history" },
  { file: "results.html",     title: "Validation record" },
  { file: "lessons.html",     title: "Bugs and lessons" },
  { file: "limits.html",      title: "Limitations" },
  { file: "roadmap.html",     title: "Status and roadmap" },
  { file: "scenario.html",    title: "Reference scenario" },
  { group: "Reference" },
  { file: "api.html",         title: "Function reference" },
  { file: "formulas.html",    title: "Formulas" },
  { file: "glossary.html",    title: "Glossary" },
  { file: "maintenance.html", title: "Docs maintenance" },
];

(function () {
  const here = location.pathname.split("/").pop() || "index.html";

  // theme: stored per viewer, optional
  function getTheme() { try { return localStorage.getItem("psrdocs-theme"); } catch (e) { return null; } }
  function setTheme(t) { try { localStorage.setItem("psrdocs-theme", t); } catch (e) {} }
  const saved = getTheme();
  if (saved) document.documentElement.setAttribute("data-theme", saved);

  // sidebar
  const side = document.createElement("nav");
  side.className = "side";
  let html = '<a class="brand" href="index.html">Pulsar navigation</a>' +
             '<div class="tagline">PulsarSimMatlab documentation</div>';
  let open = false;
  for (const p of PAGES) {
    if (p.group) {
      if (open) html += "</ul>";
      html += "<h4>" + p.group + "</h4><ul>";
      open = true;
    } else {
      const cur = p.file === here ? ' class="current"' : "";
      html += '<li><a href="' + p.file + '"' + cur + ">" + p.title + "</a></li>";
    }
  }
  if (open) html += "</ul>";
  html += '<button class="theme-toggle" type="button">Light / dark</button>';
  side.innerHTML = html;
  document.body.prepend(side);

  side.querySelector(".theme-toggle").addEventListener("click", function () {
    const cur = document.documentElement.getAttribute("data-theme") ||
      (matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light");
    const next = cur === "dark" ? "light" : "dark";
    document.documentElement.setAttribute("data-theme", next);
    setTheme(next);
  });

  // mobile menu button
  const btn = document.createElement("button");
  btn.className = "menu-btn";
  btn.type = "button";
  btn.textContent = "Menu";
  btn.addEventListener("click", function () { document.body.classList.toggle("nav-open"); });
  document.body.prepend(btn);

  // prev / next
  const pages = PAGES.filter(function (p) { return p.file; });
  const i = pages.findIndex(function (p) { return p.file === here; });
  const main = document.querySelector("main");
  if (main && i >= 0) {
    const pager = document.createElement("div");
    pager.className = "pager";
    const prev = i > 0 ? '<a href="' + pages[i - 1].file + '">← ' + pages[i - 1].title + "</a>" : "<span></span>";
    const next = i < pages.length - 1 ? '<a href="' + pages[i + 1].file + '">' + pages[i + 1].title + " →</a>" : "<span></span>";
    pager.innerHTML = prev + next;
    main.appendChild(pager);
  }

  // equations: KaTeX from vendor/katex (a local copy, so file:// still works).
  // Write \[ ... \] for a display equation and \( ... \) for inline math.
  // No $ delimiters: a dollar sign in running text must stay text.
  function load(tag, attrs, onload) {
    const el = document.createElement(tag);
    for (const k in attrs) el.setAttribute(k, attrs[k]);
    if (onload) el.addEventListener("load", onload);
    document.head.appendChild(el);
  }
  if (main && /\\\(|\\\[/.test(main.innerHTML)) {
    load("link", { rel: "stylesheet", href: "vendor/katex/katex.min.css" });
    load("script", { src: "vendor/katex/katex.min.js" }, function () {
      load("script", { src: "vendor/katex/auto-render.min.js" }, function () {
        renderMathInElement(main, {
          delimiters: [
            { left: "\\[", right: "\\]", display: true },
            { left: "\\(", right: "\\)", display: false },
          ],
          throwOnError: false,
        });
      });
    });
  }
})();
