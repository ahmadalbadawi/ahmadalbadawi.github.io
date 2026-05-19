// Theme toggle. Pairs with the inline bootstrap in _includes/theme-init.html.
// Defaults to dark. Light is opt-in. Persists across pages via two layers:
//   - localStorage["theme"]
//   - cookie "theme" (one-year max-age, SameSite=Lax)
// Either source on the next page load applies the light theme.

(function () {
  'use strict';

  var root = document.documentElement;
  var THEME_COLOR_DARK = '#10141c';
  var THEME_COLOR_LIGHT = '#fafafa';
  var COOKIE_MAX_AGE = 60 * 60 * 24 * 365; // one year

  function currentTheme() {
    return root.getAttribute('data-theme') === 'light' ? 'light' : 'dark';
  }

  function writeStorage(theme) {
    try {
      if (theme === 'light') {
        localStorage.setItem('theme', 'light');
      } else {
        localStorage.removeItem('theme');
      }
    } catch (e) {}
  }

  function writeCookie(theme) {
    if (theme === 'light') {
      document.cookie =
        'theme=light; Max-Age=' + COOKIE_MAX_AGE + '; Path=/; SameSite=Lax';
    } else {
      document.cookie = 'theme=; Max-Age=0; Path=/; SameSite=Lax';
    }
  }

  function applyTheme(theme) {
    if (theme === 'light') {
      root.setAttribute('data-theme', 'light');
    } else {
      root.removeAttribute('data-theme');
    }
    writeStorage(theme);
    writeCookie(theme);
    syncMetaThemeColor();
  }

  function syncMetaThemeColor() {
    var metas = document.querySelectorAll('meta[name="theme-color"]');
    var color = currentTheme() === 'light' ? THEME_COLOR_LIGHT : THEME_COLOR_DARK;
    metas.forEach(function (m) {
      m.removeAttribute('media');
      m.setAttribute('content', color);
    });
  }

  function syncPressed(buttons) {
    var isLight = currentTheme() === 'light';
    buttons.forEach(function (btn) {
      btn.setAttribute('aria-pressed', String(!isLight));
    });
  }

  function wire() {
    syncMetaThemeColor();
    var buttons = document.querySelectorAll('[data-theme-toggle]');
    syncPressed(buttons);

    buttons.forEach(function (btn) {
      btn.addEventListener('click', function () {
        applyTheme(currentTheme() === 'light' ? 'dark' : 'light');
        syncPressed(buttons);
      });
    });

    window.addEventListener('storage', function (e) {
      if (e.key !== 'theme') return;
      applyTheme(e.newValue === 'light' ? 'light' : 'dark');
      syncPressed(buttons);
    });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', wire);
  } else {
    wire();
  }
})();
