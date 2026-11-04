// Publications filter: by kind (journal, conference, ...) and by words in the title, authors or
// venue. Built here so the page shows no dead controls without JavaScript. The state is kept in
// the address (?type=journal&q=gpu) so a filtered list can be linked.
(function () {
  var pubs = Array.prototype.slice.call(document.querySelectorAll('.pub[data-group]'));
  if (!pubs.length) return;
  var page = document.querySelector('.page__body');
  var labels = { journal: 'Journals', conference: 'Conferences', workshop: 'Workshops', preprint: 'Preprints', other: 'Other', patent: 'Patents' };
  var counts = {};
  pubs.forEach(function (p) {
    var g = p.getAttribute('data-group');
    counts[g] = (counts[g] || 0) + 1;
    p._text = Array.prototype.map.call(p.querySelectorAll('.pub__title, .pub__authors, .pub__where'), function (e) { return e.textContent; }).join(' ').toLowerCase();
  });

  var form = document.createElement('form');
  form.className = 'pub-filter';
  form.setAttribute('role', 'search');
  form.setAttribute('aria-label', 'Filter publications');
  var opts = '<label class="pub-filter__opt"><input type="radio" name="type" value="" checked><span>All <span class="n">' + pubs.length + '</span></span></label>';
  Object.keys(labels).forEach(function (g) {
    if (!counts[g]) return;
    opts += '<label class="pub-filter__opt"><input type="radio" name="type" value="' + g + '"><span>' + labels[g] + ' <span class="n">' + counts[g] + '</span></span></label>';
  });
  form.innerHTML =
    '<fieldset><legend class="label">Show</legend>' + opts + '</fieldset>' +
    '<div class="pub-filter__search"><label class="label" for="pub-q">Search</label>' +
    '<input id="pub-q" name="q" type="search" autocomplete="off" placeholder="title, author or venue"></div>' +
    '<p class="pub-filter__count" aria-live="polite"></p>';
  form.addEventListener('submit', function (e) { e.preventDefault(); });

  var intro = page.querySelector('.intro');
  if (intro) intro.after(form); else page.prepend(form);
  var count = form.querySelector('.pub-filter__count');
  var q = form.querySelector('#pub-q');

  var params = new URLSearchParams(location.search);
  if (params.get('type') && counts[params.get('type')]) form.querySelector('input[value="' + params.get('type') + '"]').checked = true;
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
  apply();
})();
