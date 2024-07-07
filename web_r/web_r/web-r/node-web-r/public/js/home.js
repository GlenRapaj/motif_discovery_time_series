// JavaScript for carousel functionality (if needed)
// Example of dynamically loading testimonials
const testimonials = [
    { quote: "Customer quote 1.", name: "John Doe" },
    { quote: "Customer quote 2.", name: "Jane Smith" }
];

const testimonialContainer = document.getElementById('testimonial-carousel');

testimonials.forEach(testimonial => {
    const div = document.createElement('div');
    div.classList.add('testimonial');
    div.innerHTML = `
        <p>"${testimonial.quote}"</p>
        <cite>- ${testimonial.name}</cite>
    `;
    testimonialContainer.appendChild(div);
});