// Theme toggle - pairs with the inline bootstrap in _includes/theme-init.html.
// The site defaults to dark; clicking the button switches between dark and
// light and persists "light" in localStorage. Switching back to dark removes
// the saved value so the default takes over.

(function () {
  const root = document.documentElement;
  const THEME_COLOR_DARK = '#10141c';
  const THEME_COLOR_LIGHT = '#fafafa';

  function currentTheme() {
    return root.getAttribute('data-theme') === 'light' ? 'light' : 'dark';
  }

  function applyTheme(theme) {
    if (theme === 'light') {
      root.setAttribute('data-theme', 'light');
      try { localStorage.setItem('theme', 'light'); } catch (e) {}
    } else {
      root.removeAttribute('data-theme');
      try { localStorage.removeItem('theme'); } catch (e) {}
    }
    syncMetaThemeColor();
  }

  function syncMetaThemeColor() {
    const metas = document.querySelectorAll('meta[name="theme-color"]');
    const color = currentTheme() === 'light' ? THEME_COLOR_LIGHT : THEME_COLOR_DARK;
    metas.forEach((m) => {
      // Drop any media-query qualifier so the value applies unconditionally.
      m.removeAttribute('media');
      m.setAttribute('content', color);
    });
  }

  function syncPressed(buttons) {
    const isLight = currentTheme() === 'light';
    buttons.forEach((btn) => btn.setAttribute('aria-pressed', String(!isLight)));
  }

  document.addEventListener('DOMContentLoaded', () => {
    syncMetaThemeColor();
    const buttons = document.querySelectorAll('[data-theme-toggle]');
    syncPressed(buttons);

    buttons.forEach((btn) => {
      btn.addEventListener('click', () => {
        applyTheme(currentTheme() === 'light' ? 'dark' : 'light');
        syncPressed(buttons);
      });
    });

    // Cross-tab sync.
    window.addEventListener('storage', (e) => {
      if (e.key !== 'theme') return;
      applyTheme(e.newValue === 'light' ? 'light' : 'dark');
      syncPressed(buttons);
    });
  });
})();
