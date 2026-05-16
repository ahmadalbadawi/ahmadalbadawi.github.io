// Mobile nav toggle - no jQuery. Falls back to CSS-only behavior if JS fails.

(function () {
  document.addEventListener('DOMContentLoaded', () => {
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

    // Close on Escape.
    document.addEventListener('keydown', (e) => {
      if (e.key === 'Escape' && nav.dataset.open === 'true') {
        setOpen(false);
        toggle.focus();
      }
    });

    // Close on outside click.
    document.addEventListener('click', (e) => {
      if (nav.dataset.open !== 'true') return;
      if (nav.contains(e.target)) return;
      setOpen(false);
    });

    // Close when crossing the desktop breakpoint.
    const mql = window.matchMedia('(min-width: 768px)');
    mql.addEventListener('change', (e) => {
      if (e.matches) setOpen(false);
    });
  });
})();
