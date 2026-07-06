Interactive diagrams
====================

Two self-contained, interactive diagrams of how TOPO fits together. Each is a
full-screen SVG with clickable nodes, flow chips that light up and animate a
sequence, and a light/dark theme toggle. They open in a new tab.

.. rubric:: Architecture — the stack

The four runners (``topo-mdrun``, ``topo-optimize``, ``topo-csp``,
``topo-cylinder``), the one shared MD engine they all delegate to, the
coarse-grained force field, and the outputs. Flow chips trace *Equilibrium*,
*Optimize*, *Synthesis*, *Cylinder*, and *Analysis* paths through the stack.

`Open the architecture diagram → <architecture.html>`_

.. rubric:: Call flow — end to end

The end-to-end **function-call sequence** for the two workflows side by side:
isolated-protein MD (``topo-mdrun``) and protein synthesis
(``topo-csp``), both funneling into the shared per-segment core
(``setup_simulation`` → ``attach_reporters`` → ``step`` → ``finalize_simulation``).
Numbered badges show call order; chips replay *Isolated MD*, *CSP synthesis*, the
*per-residue* 3-stage loop, *Model build*, and *Engine step*.

`Open the call-flow diagram → <callflow.html>`_

.. raw:: html

   <style>
     .diagram-cards { display: flex; flex-wrap: wrap; gap: 16px; margin: 1.5rem 0 0.5rem; }
     .diagram-cards a.dcard {
       flex: 1 1 260px; display: block; padding: 16px 18px; border-radius: 12px;
       border: 1px solid var(--color-background-border, #ccc);
       background: var(--color-background-secondary, #f5f5f5);
       text-decoration: none; transition: border-color .15s ease, transform .15s ease;
     }
     .diagram-cards a.dcard:hover { border-color: var(--color-brand-primary, #d97757); transform: translateY(-2px); }
     .diagram-cards .dcard .dtitle { font-weight: 700; font-size: 1.05rem; margin-bottom: 4px;
       color: var(--color-brand-content, var(--color-content-foreground)); }
     .diagram-cards .dcard .ddesc { font-size: 0.86rem; color: var(--color-foreground-secondary, #666); }
   </style>
   <div class="diagram-cards">
     <a class="dcard" href="architecture.html" target="_blank" rel="noopener">
       <div class="dtitle">Architecture — the stack &#8599;</div>
       <div class="ddesc">Four runners · one MD engine · one force field. Animated Equilibrium / Optimize / Synthesis / Cylinder / Analysis paths.</div>
     </a>
     <a class="dcard" href="callflow.html" target="_blank" rel="noopener">
       <div class="dtitle">Call flow — end to end &#8599;</div>
       <div class="ddesc">Isolated-protein MD and CSP synthesis as an ordered function-call sequence, sharing the engine core.</div>
     </a>
   </div>
