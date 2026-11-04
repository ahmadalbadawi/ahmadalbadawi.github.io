// Publications filter: by kind (journal, conference, ...) and by words in the title, authors or
// venue. The controls are rendered by _includes/pub-filter.html; this only wires them up. The
// state is kept in the address (?type=journal&q=gpu) so a filtered list can be linked.
(function () {
  var form = document.querySelector('.pub-filter');
  var pubs = Array.prototype.slice.call(document.querySelectorAll('.pub[data-group]'));
  if (!form || !pubs.length) return;
  var page = document.querySelector('.page__body');
  var count = form.querySelector('.pub-filter__count');
  var q = form.querySelector('#pub-q');
  pubs.forEach(function (p) {
    p._text = Array.prototype.map.call(p.querySelectorAll('.pub__title, .pub__authors, .pub__where'), function (e) { return e.textContent; }).join(' ').toLowerCase();
  });
  form.addEventListener('submit', function (e) { e.preventDefault(); });

  var params = new URLSearchParams(location.search);
  var preset = params.get('type') && form.querySelector('input[name="type"][value="' + params.get('type').replace(/[^a-z]/g, '') + '"]');
  if (preset) preset.checked = true;
  if (params.get('q')) q.value = params.get('q');

  function apply() {
    var type = form.querySelector('input[name="type"]:checked').value;
    var words = q.value.toLowerCase().split(/\s+/).filter(Boolean);
    var shown = 0;
    pubs.forEach(function (p) {
      var ok = (!type || p.getAttribute('data-group') === type) &&
        words.every(function (w) { return p._text.indexOf(w) !== -1; });
      p.hidden = !ok;
      if (ok) shown++;
    });
    // Hide a year (or the patents heading) once nothing under it is left.
    page.querySelectorAll('.pub-list').forEach(function (list) {
      var empty = !list.querySelector('.pub:not([hidden])');
      var section = list.closest('.year-section');
      if (section) { section.hidden = empty; return; }
      list.hidden = empty;
      var h = list.previousElementSibling;
      if (h && h.tagName === 'H2') h.hidden = empty;
    });
    count.textContent = shown === pubs.length ? pubs.length + ' entries' : 'Showing ' + shown + ' of ' + pubs.length;

    var next = new URLSearchParams();
    if (type) next.set('type', type);
    if (q.value.trim()) next.set('q', q.value.trim());
    var s = next.toString();
    history.replaceState(null, '', location.pathname + (s ? '?' + s : ''));
  }
  form.addEventListener('change', apply);
  q.addEventListener('input', apply);
  if (preset || q.value) apply();
})();
