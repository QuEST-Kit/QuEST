
window.addEventListener('DOMContentLoaded', (event) => {

    /* injects the QuEST version number below the logo in the 
    * left nav-bar; this requires this script is included in
    * header.html, and the version number style is customised
    * in *-quest.css. This snippet was taken from Zephyr's
    * github (github.com/zephyrproject-rtos/zephyr) at:
    * zephyr/doc/_doxygen/doxygen-awesome-zephyr.js */

    let version = document.getElementById('projectnumber').innerText
    let titleTable = document.querySelector('#titlearea table');
    let cell = titleTable.insertRow(1).insertCell(0);
    cell.innerHTML = '<div id="projectversion">' + version + '</div>';


    /* deletes the superfluous/ugly title atop Mainpage (which
    * says "The Quantum Exact Simulation Toolkit"), before the
    * banner image. The same html element is necessary on other
    * pages, so we selectively delete it from index.html. This
    * would be more elegantly done using css in *-quest.css, 
    * which would avoid the post-render-delete visual stutter,
    * but alas I cannot find unique page-ids I can refer to! */

    let fn = window.location.pathname.split("/").pop();
    if (fn === "index.html") {
        const element = document.querySelector(".header");
        if (element) {
            element.remove();
        }
    }

});
