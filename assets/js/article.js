// Posts only. Everything here is an enhancement: without it the contents list is a disclosure,
// footnotes are end notes and code blocks scroll.
(function () {
  var wideToc = matchMedia('(min-width: 80rem)');
  var wideMargin = matchMedia('(min-width: 72rem)');
  var body = document.getElementById('post-body');
  if (!body) return;

  // Contents: a fixed list beside the text on wide screens, a closed disclosure above it
  // otherwise. Marks the section being read.
  var toc = document.querySelector('.toc');
  if (toc) {
    var details = toc.querySelector('details');
    toc.setAttribute('data-js', '');
    var setOpen = function () { details.open = wideToc.matches; };
    setOpen();
    wideToc.addEventListener('change', setOpen);

    var links = {};
    var targets = [];
    toc.querySelectorAll('details > ul > li > a[href^="#"]').forEach(function (a) {
      var el = document.getElementById(decodeURIComponent(a.getAttribute('href').slice(1)));
      if (el) { links[el.id] = a; targets.push(el); }
    });
    if ('IntersectionObserver' in window && targets.length) {
      var io = new IntersectionObserver(function (entries) {
        entries.forEach(function (e) {
          if (!e.isIntersecting) return;
          toc.querySelectorAll('[aria-current]').forEach(function (a) { a.removeAttribute('aria-current'); });
          links[e.target.id].setAttribute('aria-current', 'location');
        });
      }, { rootMargin: '0px 0px -70% 0px' });
      targets.forEach(function (t) { io.observe(t); });
    }
  }

  // Sidenotes: on wide screens each kramdown footnote moves beside its reference, so a reader
  // meets it once, where it is cited. On narrow screens it goes back to the end notes.
  var notes = body.querySelector('.footnotes');
  var moved = [];
  function toMargin() {
    if (moved.length) return;
    body.querySelectorAll('a.footnote[href^="#fn"]').forEach(function (a) {
      var li = document.getElementById(decodeURIComponent(a.getAttribute('href').slice(1)));
      var sup = a.closest('sup');
      if (!li || !sup) return;
      var n = a.textContent.trim();
      var note = document.createElement('span');
      note.className = 'sidenote';
      note.setAttribute('role', 'note');
      note.innerHTML = '<span class="visually-hidden">Note ' + n + ': </span><span class="n" aria-hidden="true">' + n + '</span>';
      while (li.firstChild) note.appendChild(li.firstChild);
      note.querySelectorAll('.reversefootnote').forEach(function (b) { b.hidden = true; });
      sup.classList.add('sn-ref');
      a.setAttribute('aria-hidden', 'true');
      a.tabIndex = -1;
      sup.after(note);
      moved.push({ a: a, li: li, note: note });
    });
    if (moved.length && notes) notes.hidden = true;
  }
  function toEnd() {
    moved.forEach(function (m) {
      m.note.querySelectorAll('.reversefootnote').forEach(function (b) { b.hidden = false; });
      m.note.firstChild.remove();
      m.note.firstChild.remove();
      while (m.note.firstChild) m.li.appendChild(m.note.firstChild);
      m.note.remove();
      m.a.removeAttribute('aria-hidden');
      m.a.removeAttribute('tabindex');
      m.a.closest('sup').classList.remove('sn-ref');
    });
    moved = [];
    if (notes) notes.hidden = false;
  }
  if (notes) {
    if (wideMargin.matches) toMargin();
    wideMargin.addEventListener('change', function (e) { if (e.matches) toMargin(); else toEnd(); });
  }

  // Code: scrollable blocks are reachable by keyboard, and each has a Copy button.
  body.querySelectorAll('pre').forEach(function (pre) { pre.tabIndex = 0; });
  focusableTables();
  if (navigator.clipboard) {
    body.querySelectorAll('div.highlighter-rouge').forEach(function (block) {
      var code = block.querySelector('code');
      if (!code) return;
      var btn = document.createElement('button');
      btn.type = 'button';
      btn.className = 'copy-code';
      btn.textContent = 'Copy';
      btn.setAttribute('aria-label', 'Copy code to clipboard');
      btn.addEventListener('click', function () {
        navigator.clipboard.writeText(code.innerText.replace(/\n$/, '')).then(function () {
          btn.textContent = 'Copied';
          setTimeout(function () { btn.textContent = 'Copy'; }, 2000);
        });
      });
      block.appendChild(btn);
    });
  }

  // Tables wider than the column scroll (CSS); make those reachable by keyboard.
  function focusableTables() {
    body.querySelectorAll('table').forEach(function (t) { if (t.scrollWidth > t.clientWidth + 1) t.tabIndex = 0; });
  }

  // Maths wider than the column scrolls inside its own box instead of widening the page, and
  // can then be reached by keyboard. Display maths scrolls by CSS; a rare inline formula wider
  // than the column (only on phones) is turned into a scrolling block here.
  function fitMath() {
    body.querySelectorAll('mjx-container:not([display="true"])').forEach(function (m) {
      var block = m.parentElement.closest('p, li, td, th, blockquote, div');
      if (block && m.getBoundingClientRect().width > block.clientWidth + 1) m.classList.add('mjx-scroll');
    });
    body.querySelectorAll('mjx-container[display="true"], mjx-container.mjx-scroll').forEach(function (m) {
      if (m.scrollWidth > m.clientWidth + 1) m.tabIndex = 0;
    });
    focusableTables();
  }
  var waitForMathJax = setInterval(function () {
    if (window.MathJax && MathJax.startup && MathJax.startup.promise) {
      clearInterval(waitForMathJax);
      MathJax.startup.promise.then(fitMath);
    }
  }, 200);
  setTimeout(function () { clearInterval(waitForMathJax); }, 30000);
})();
