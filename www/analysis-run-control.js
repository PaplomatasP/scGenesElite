(function () {
  "use strict";

  // The R process is blocked for the whole analysis, so the swap to
  // "Stop analysis" is done here instead of waiting for a server round trip.
  var RUNNING = "analysis-running";
  var STOPPING = "analysis-stopping";
  var BUSY = "app-busy";
  var BUSY_DELAY_MS = 400;
  var busyTimer = null;

  function setState(running, stopping) {
    var body = document.body;
    if (!body) return;
    body.classList.toggle(RUNNING, !!running);
    body.classList.toggle(STOPPING, !!running && !!stopping);
  }

  function onRunClicked() {
    setState(true, false);
  }

  // The click is queued until R is free again, so lock the button meanwhile:
  // a stop cannot be fired repeatedly either.
  function onStopClicked() {
    setState(true, true);
  }

  document.addEventListener("click", function (event) {
    var target = event.target;
    if (!target || !target.closest) return;
    if (target.closest("#click")) return onRunClicked();
    if (target.closest("#stopAnalysis")) return onStopClicked();
  });

  // A slim sticky banner shown for as long as R is occupied, injected once the
  // page body exists. Every plot/table re-render counts, not just the initial
  // run, since changing a setting (e.g. turning the heatmap on) recomputes it.
  function ensureBusyBar() {
    if (document.querySelector(".app-busy-bar")) return;
    var bar = document.createElement("div");
    bar.className = "app-busy-bar";
    bar.innerHTML = '<span class="spinner" aria-hidden="true"></span>' +
      '<span>Still working - this can take a while for larger datasets or heavier methods.</span>';
    document.body.insertBefore(bar, document.body.firstChild);
  }

  // A fullscreen view for any result card: add class="fullscreen-btn" and
  // data-target="<the plotOutput/tableOutput id>" to a button, it clones
  // whatever is currently rendered there into a large centered overlay.
  function ensureLightbox() {
    if (document.querySelector(".sc-lightbox")) return;
    var box = document.createElement("div");
    box.className = "sc-lightbox";
    box.innerHTML =
      '<div class="sc-lightbox-panel">' +
      '<button type="button" class="sc-lightbox-close" aria-label="Close">&times;</button>' +
      '<div class="sc-lightbox-body"></div>' +
      "</div>";
    document.body.appendChild(box);

    function close() {
      box.classList.remove("is-open");
      box.querySelector(".sc-lightbox-body").innerHTML = "";
    }
    box.addEventListener("click", function (event) {
      if (event.target === box || event.target.closest(".sc-lightbox-close")) close();
    });
    document.addEventListener("keydown", function (event) {
      if (event.key === "Escape") close();
    });
  }

  function openFullscreen(targetId) {
    var source = document.getElementById(targetId);
    var box = document.querySelector(".sc-lightbox");
    if (!source || !box) return;
    var body = box.querySelector(".sc-lightbox-body");
    body.innerHTML = source.innerHTML;
    box.classList.add("is-open");
  }

  document.addEventListener("click", function (event) {
    var trigger = event.target.closest && event.target.closest(".fullscreen-btn");
    if (trigger && trigger.dataset.target) openFullscreen(trigger.dataset.target);
  });

  var $ = window.jQuery;
  if ($) {
    // Fires whenever Shiny sends the button value, whatever triggered it.
    $(document).on("shiny:inputchanged", function (event) {
      if (event.name === "click") onRunClicked();
      else if (event.name === "stopAnalysis") onStopClicked();
    });

    // "Please wait" has to outlive MethodData(): the plots, heatmaps and
    // networks render after its withProgress() has already closed, and the app
    // is frozen for all of it. shiny:busy/idle span the whole round trip, so
    // this also covers a slow re-render after a settings change.
    $(document).on("shiny:busy", function () {
      if (busyTimer !== null) return;
      busyTimer = window.setTimeout(function () {
        document.body.classList.add(BUSY);
      }, BUSY_DELAY_MS);
    });
    $(document).on("shiny:idle", function () {
      window.clearTimeout(busyTimer);
      busyTimer = null;
      document.body.classList.remove(BUSY);
    });
  }

  // The server resets the pair once the run has released the process.
  function registerHandler() {
    if (!window.Shiny || !Shiny.addCustomMessageHandler) return false;
    Shiny.addCustomMessageHandler("scgenes-analysis-state", function (message) {
      setState(message && message.running, message && message.stopping);
    });
    return true;
  }

  function init() {
    ensureBusyBar();
    ensureLightbox();
    if (!registerHandler()) document.addEventListener("DOMContentLoaded", registerHandler);
  }

  if (document.body) init();
  else document.addEventListener("DOMContentLoaded", init);
})();
