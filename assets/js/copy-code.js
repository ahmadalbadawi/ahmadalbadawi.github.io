// Adds a "Copy" button to each Rouge-highlighted code block.

(function () {
  function init() {
    const blocks = document.querySelectorAll('.highlight');
    blocks.forEach((block) => {
      // Wrap so the button can position relative to the block.
      const wrap = document.createElement('div');
      wrap.className = 'code-block-wrap';
      block.parentNode.insertBefore(wrap, block);
      wrap.appendChild(block);

      const btn = document.createElement('button');
      btn.type = 'button';
      btn.className = 'copy-code-btn';
      btn.setAttribute('aria-label', 'Copy code to clipboard');
      btn.textContent = 'Copy';
      wrap.appendChild(btn);

      btn.addEventListener('click', () => {
        const code = block.querySelector('code');
        if (!code) return;
        navigator.clipboard.writeText(code.innerText.trim()).then(() => {
          btn.textContent = 'Copied';
          setTimeout(() => { btn.textContent = 'Copy'; }, 2500);
        });
      });
    });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
