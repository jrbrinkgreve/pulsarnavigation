// Pulsar navigation docs: sidebar, prev/next links and theme toggle.
// The page list lives here only, so adding a page = one line in PAGES.
// Works from file:// (no server, no external resources).

const PAGES = [
  { group: "Start" },
  { file: "index.html",       title: "Home" },
  { file: "physics.html",     title: "Physics background" },
  { file: "conventions.html", title: "Conventions and why" },
  { file: "pipeline.html",    title: "Pipeline overview" },
  { file: "parameters.html",  title: "Parameters" },
  { group: "Stages" },
  { file: "sky.html",         title: "1–2 Sky: pulsar and dispersion" },
  { file: "receiver.html",    title: "3–4 Receiver: noise, RFI, IQ" },
  { file: "fullband.html",    title: "5–6 Full-band front end" },
  { file: "channels.html",    title: "5b–6b Channelized front end" },
  { file: "excision.html",    title: "RFI excision" },
  { file: "fold.html",        title: "7 Folding" },
  { file: "toa.html",         title: "8 TOA estimation" },
  { file: "detection.html",   title: "Detection and good TOAs" },
  { file: "validation.html",  title: "9 Validation and checks" },
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
})();
