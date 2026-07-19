// Collapse each release entry on the Changelog page into a native <details> so a
// long history (10+ releases) stays scannable — you see a list of release headings
// and expand the one you want. Runs ONLY on the changelog page; CHANGELOG.md stays
// clean Markdown (single source of truth) — the collapsing is purely presentational.
document.addEventListener("DOMContentLoaded", function () {
  var root = document.getElementById("changelog");
  if (!root) return; // not the changelog page

  // Each release is a level-2 <section> directly under the H1 "Changelog" section;
  // its ### subsections (Added / Changed / …) nest inside it, so wrapping the whole
  // section collapses the entire entry.
  root.querySelectorAll(":scope > section").forEach(function (sec) {
    var h2 = sec.querySelector(":scope > h2");
    if (!h2) return;

    var details = document.createElement("details");
    details.className = "changelog-entry";

    var summary = document.createElement("summary");
    while (h2.firstChild) summary.appendChild(h2.firstChild); // heading text + ¶ link
    h2.remove();
    details.appendChild(summary);

    while (sec.firstChild) details.appendChild(sec.firstChild); // the rest of the entry
    sec.appendChild(details);
  });

  // If the page was opened with a fragment pointing inside a collapsed entry
  // (e.g. a deep link to #id2 or #added), open that entry and scroll to it.
  function openHashTarget() {
    if (!location.hash) return;
    var target = document.getElementById(decodeURIComponent(location.hash.slice(1)));
    if (!target) return;
    var d = target.closest("details.changelog-entry");
    if (d && !d.open) {
      d.open = true;
      target.scrollIntoView();
    }
  }
  openHashTarget();
  window.addEventListener("hashchange", openHashTarget);
});
