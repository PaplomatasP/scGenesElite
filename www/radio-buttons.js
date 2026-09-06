(function () {
  "use strict";

  function closeMobileNavigation() {
    var collapse = document.querySelector(".navbar-collapse");
    if (!collapse || !window.jQuery) return;

    var $collapse = window.jQuery(collapse);
    if ($collapse.hasClass("in")) {
      $collapse.collapse("hide");
    } else if ($collapse.hasClass("collapsing")) {
      $collapse.one("shown.bs.collapse", function () {
        $collapse.collapse("hide");
      });
    }
  }

  function navigateToTab(value) {
    if (!window.jQuery || !value) return;

    var target = window.jQuery('.navbar-nav a[data-value="' + value + '"]').first();
    if (!target.length) return;

    target.tab("show");
    closeMobileNavigation();
    window.requestAnimationFrame(function () {
      window.scrollTo({ top: 0, behavior: "smooth" });
    });
  }

  document.addEventListener("click", function (event) {
    var navTrigger = event.target.closest(".sc-nav-link");
    if (navTrigger) {
      event.preventDefault();
      navigateToTab(navTrigger.getAttribute("data-nav-target"));
      return;
    }

    var navLink = event.target.closest(".navbar-nav a");
    if (navLink && !navLink.classList.contains("dropdown-toggle")) {
      closeMobileNavigation();
    }
  });

  if (window.jQuery) {
    window.jQuery(document).on("shown.bs.tab", ".navbar-nav a[data-toggle='tab']", function () {
      var item = window.jQuery(this).closest("li");
      var parentDropdown = item.closest("li.dropdown");
      window.jQuery(".navbar-nav > li.dropdown").removeClass("active");
      if (parentDropdown.length) parentDropdown.addClass("active");
    });
  }
})();
