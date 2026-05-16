// Mobile nav toggle + desktop "More" dropdowns. No jQuery. Falls back to
// CSS-only behavior if JS fails (see noscript styles in _includes/head.html).

(function () {
  document.addEventListener('DOMContentLoaded', () => {
    initMobileToggle();
    initDropdowns();
  });

  function initMobileToggle() {
    const nav = document.querySelector('.site-nav');
    const toggle = document.querySelector('[data-nav-toggle]');
    if (!nav || !toggle) return;

    function setOpen(open) {
      nav.dataset.open = open ? 'true' : 'false';
      toggle.setAttribute('aria-expanded', String(open));
    }

    toggle.addEventListener('click', () => {
      setOpen(nav.dataset.open !== 'true');
    });

    document.addEventListener('keydown', (e) => {
      if (e.key === 'Escape' && nav.dataset.open === 'true') {
        setOpen(false);
        toggle.focus();
      }
    });

    document.addEventListener('click', (e) => {
      if (nav.dataset.open !== 'true') return;
      if (nav.contains(e.target)) return;
      setOpen(false);
    });

    const mql = window.matchMedia('(min-width: 768px)');
    mql.addEventListener('change', (e) => {
      if (e.matches) setOpen(false);
    });
  }

  function initDropdowns() {
    const dropdowns = document.querySelectorAll('[data-nav-dropdown]');
    if (!dropdowns.length) return;
    const desktop = window.matchMedia('(min-width: 768px)');

    function closeAll(exceptEl) {
      dropdowns.forEach((d) => {
        if (d === exceptEl) return;
        const btn = d.querySelector('button');
        const menu = d.querySelector('ul');
        if (btn) btn.setAttribute('aria-expanded', 'false');
        if (menu) menu.hidden = true;
      });
    }

    dropdowns.forEach((d) => {
      const btn = d.querySelector('button');
      const menu = d.querySelector('ul');
      if (!btn || !menu) return;

      btn.addEventListener('click', (e) => {
        if (!desktop.matches) return; // mobile flattens the menu, no toggle
        e.stopPropagation();
        const isOpen = btn.getAttribute('aria-expanded') === 'true';
        closeAll(isOpen ? null : d);
        btn.setAttribute('aria-expanded', String(!isOpen));
        menu.hidden = isOpen;
      });
    });

    document.addEventListener('click', () => closeAll());

    document.addEventListener('keydown', (e) => {
      if (e.key === 'Escape') closeAll();
    });

    desktop.addEventListener('change', () => closeAll());
  }
})();
