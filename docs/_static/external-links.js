// Open external links (different host — e.g. DOI links, GitHub, OpenMM) in a new
// browser tab. Sphinx/reST renders bare URLs as plain <a> tags with no target,
// so this runs after the page loads and adds target="_blank" (plus the
// rel="noopener noreferrer" security hardening) to every off-site link.
document.addEventListener("DOMContentLoaded", function () {
  var here = window.location.hostname;
  document.querySelectorAll("a[href^='http']").forEach(function (a) {
    if (a.hostname && a.hostname !== here) {
      a.setAttribute("target", "_blank");
      a.setAttribute("rel", "noopener noreferrer");
    }
  });
});
